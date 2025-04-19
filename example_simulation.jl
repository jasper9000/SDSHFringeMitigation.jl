### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ d09ffcda-1c5e-11f0-2a95-13c00027ecfa
begin
	# push!(LOAD_PATH, dirname(pwd()));
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ 788fdd18-cf33-40f5-8845-9a35bed9c34d
using Revise


# ╔═╡ 89282ccd-b3c4-4587-af1f-96b6582b820a
using SDSHFringeMitigation; const SFM = SDSHFringeMitigation;

# ╔═╡ 111f1904-5d2d-41ea-ad40-8a1b85dcf37b
begin
	import Plots
	Plots.default(dpi=300)
	Plots.pythonplot()
	palette = :Set1_6
	# import PlotlyJS as pltjs
	# Plots.plotlyjs()
	# using Plots
	# default(dpi=300)
	# plotlyjs()
	
	# using PlutoPlotly
	using LaTeXStrings, ColorSchemes, ProgressBars, PlutoUI, Suppressor
	# import PlutoUIExtra
	using PartialFunctions, Statistics, LinearAlgebra, DSP
end

# ╔═╡ c27eafaf-cb46-4e90-9492-945221ebcca7
include("src/SDSHFringeMitigation.jl");

# ╔═╡ 0742953a-7f12-4a94-a42a-9433544bcfc1
HTML("""
<!-- the wrapper span -->
<div>
	<button id="myrestart" href="#">Restart Notebook</button>
	
	<script>
		const div = currentScript.parentElement
		const button = div.querySelector("button#myrestart")
		console.log(button);
		button.onclick = function() { restart_nb() };
		console.log(div.dispatchEvent)
		function restart_nb() {
			console.log("Send event");
			editor_state.notebook.process_status = "no_process";
			window.dispatchEvent(
				new CustomEvent("restart_process"),
				{},
				{ notebook_id: editor_state.notebook.id }
			);
		};
	</script>
</div>
""")

# ╔═╡ a6eee18b-7672-4154-88c8-fddb2dc5b6cf
# PlutoUIExtra.Sidebar(
# 	"abc "^10,
# 	md"---",
# 	(@bind a Slider(1:10)),
# )

# ╔═╡ 5107465f-5c73-41dd-af15-df87dc0b3d76
md"""
# Simulation Signal Parameters

f\_sampling: $(@bind f_sampling_MHz NumberField(10:10:1000, default=20)) MHz

Optical Delay Path: $(@bind l_delay_m NumberField(100:10:5_000, default=500)) m
"""

# ╔═╡ df222550-e9ca-4bc0-bb01-67b6079c16ef
@bind pn_parameters PlutoUI.combine() do child
	md"""
	# Phase Noise Shape Parameters
	- Beat Note SNR: $(@bind beat_note_snr_dB Slider(0:1:50, default=30, show_value=true)) dB
	- Lorentzian Linewidth: $(child(Slider(1:0.1:5; default=2.2, show_value=true)))
	- 1/f Noise Scale: $(child(Slider(1:0.1:15; default=7.7, show_value=true)))
	- Gaussian Bump Amplitude: $(child(Slider(3:0.1:15; default=10.0, show_value=true)))
	- Gaussian Bump Center Frequency: $(child(Slider(3:0.1:7; default=6.0, show_value=true)))
	- Gaussian Bump Width: $(child(Slider(2:0.1:7; default=5.5, show_value=true)))
	"""
end

# ╔═╡ 5c394cd1-c50f-4048-ac08-c6d0e90b1851
md"""
# KRR Parameters

Kernel Function: $(@bind kernel_func Select([SFM.rbf!, SFM.matern12!, SFM.matern32!, SFM.matern52!, SFM.matern72!, SFM.matern92!]))

SNR Threshold: $(@bind min_snr_dB NumberField(1.0:1.0:30.0, default=3)) dB

N training samples: $(@bind krr_train_samples NumberField(300:100:1000, default=600))

KRR x-axis log scaling: $(@bind xlog CheckBox(default=false))
KRR y-axis log scaling: $(@bind ylog CheckBox(default=true))

Estimator from which training samples are selected: $(@bind training_estimator Select(["S_INV - P", "S_INV"]))

KRR Hyperparameter Gridsearch size: $(@bind n_krr_grid NumberField(5:5:50, default=15))
"""

# ╔═╡ 369fa5fe-f2c1-4676-ab68-c3f6a8292656
begin
	f_sampling = f_sampling_MHz * 1e6
	n_supersample = 3
	f_sampling_super = f_sampling * n_supersample
	L = 4_000_000 * n_supersample

	f_beat = f_sampling / 4
	
	# define delay
	n_si = 1.5
	c = 3e8
	tau_d = l_delay_m * n_si / c
	tau_samples = round(Int, tau_d * f_sampling_super)

	gaussian_pdf = (f, f_center, width) -> exp.(-0.5 * ((f .- f_center) ./ width).^2) ./ (width * sqrt(2 * π))

	noise_func_lorentz = (f) -> 10.0.^pn_parameters[1]./π
	noise_func_flicker = (f) -> 10.0.^pn_parameters[2] .* f.^(-1)
	noise_bump = (f) -> 10.0.^pn_parameters[3] .* gaussian_pdf(f, 10.0.^pn_parameters[4], 10.0.^pn_parameters[5])
	
	noise_func_fn = (f) -> noise_func_flicker(f) .+ noise_bump(f) .+ noise_func_lorentz(f)
	
	noise_func_pn = (f) -> noise_func_fn(f) ./ f.^2

end;

# ╔═╡ a09fa9b8-8350-4951-a133-e50f7182db79
begin
	Plots.gr()
	# ff, pp = SFM.welch(pn_custom, f_sampling_super, 10)
	# pp .*= ff.^2
	f_noise = 10.0.^ (1:0.01:7)
	
	p_pn = Plots.plot(xscale=:log10, yscale=:log10, palette=palette, leg=:bottomleft)
	Plots.plot!(p_pn, f_noise, noise_func_lorentz(f_noise) ./ f_noise.^2, label="Lorentz", linewidth=2.5, linestyle=:dash)
	Plots.plot!(p_pn, f_noise, noise_func_flicker(f_noise) ./ f_noise.^2, label="1/f Noise", linewidth=2.5, linestyle=:dash)
	Plots.plot!(p_pn, f_noise, noise_bump(f_noise) ./ f_noise.^2, label="Bump", linewidth=2.5, linestyle=:dash)
	Plots.plot!(p_pn, SFM.ds(f_noise, noise_func_pn(f_noise), xlog=true), label="Total", linewidth=2.5, ylabel=L"PN-PSD / $rad^2 Hz^{-1}$")
	Plots.hline!(p_pn, [2/(f_sampling * 10.0^(beat_note_snr_dB/10))], linewidth=2, label="Measurement Noise Floor")
	Plots.ylims!(p_pn, 1e-20, 1e5)

	
	p_fn = Plots.plot(xscale=:log10, yscale=:log10, label="True FN-PSD", palette=palette, xlabel="Frequency / Hz", ylabel=L"FN-PSD / $Hz^2 Hz^{-1}$",
					 leg=:bottomleft)
	Plots.plot!(p_fn, f_noise, noise_func_lorentz(f_noise) ./ f_noise.^0, label="Lorentz", linewidth=2.5, linestyle=:dash)
	Plots.plot!(p_fn, f_noise, noise_func_flicker(f_noise), label="1/f Noise", linewidth=2.5, linestyle=:dash)
	Plots.plot!(p_fn, f_noise, noise_bump(f_noise), label="Bump", linewidth=2.5, linestyle=:dash)
	Plots.plot!(p_fn, SFM.ds(f_noise, noise_func_fn(f_noise), xlog=true), label="Total", linewidth=2.5)
	Plots.plot!(p_fn, f_noise, f_noise.^2 .* 2/(f_sampling * 10.0^(beat_note_snr_dB/10)), linewidth=2, label="Measurement Noise Floor")
	# Plots.ylims!(1, minimum([minimum(noise_func_fn(f_noise)), 1e10]))
	Plots.ylims!(1, 1e6)
	
	Plots.plot(p_pn, p_fn, layout=(2,1), sharex=true)

	
	# Plots.plot!(p1, xlabel="Frequency / Hz", ylabel=L"PN-PSD / $rad^2\;Hz^{-1}$",
	# 			xticks=[10^i for i in 1:7],
	# 		    yticks=[10.0^i for i in -14:2:0],
	# 		   title="Simulated Phase Noise")
	# PlutoPlot(p1)
	# p1
end

# ╔═╡ 8da77ef1-ba12-4883-9f25-88fad120772a
pn_parameters; beat_note_snr_dB; f_sampling; l_delay_m; md"""
# Run everything?

The checkmark resets if phase noise parameters are changed, to save on computation.

GO! $(@bind run_all CheckBox(default=false))

"""

# ╔═╡ 81ccb6d3-6608-49d1-a554-04a11a7839b0


# ╔═╡ 2b4666d8-16cb-46d0-bd48-89450f3e70ba
begin
	if run_all
		
		pn_custom = SFM.generate_custom_noise(L, f_sampling_super, noise_func_pn);
		
		# calculate Δϕ
		pn_delayed_supersampled = pn_custom[1:end-tau_samples] .- pn_custom[tau_samples+1:end]
		# downsample
		sos = digitalfilter(Lowpass((f_sampling/f_sampling_super)/2), Butterworth(8))
		pn_delayed = filtfilt(sos, pn_delayed_supersampled)[1:n_supersample:end]
	
		# generate measurement noise
	
	
		var_mn = 3e-9
		mn_noise_function = (f) -> var_mn * 2 / f_sampling
		mn = SFM.generate_custom_noise(L, f_sampling, mn_noise_function);
	
		# define signal parameters
		beat_note_snr_lin = 10.0 ^ (beat_note_snr_dB/10.0)
		A_0 = sqrt(2 * var_mn * beat_note_snr_lin)
		
		t = collect(1:length(pn_delayed)) ./ f_sampling
		y = A_0 .* sin.(2π * f_beat * t .+ pn_delayed) .+ mn[1:length(pn_delayed)];
	end
end

# ╔═╡ 646204c0-4bde-4f69-bfa9-3c481e03fb76
begin
	if run_all
		# calculate and plot spectrum
		ff_y, pp_y = SFM.welch(y, f_sampling, 10)
		p2 = Plots.plot(palette=palette)
		Plots.plot!(p2, SFM.ds((ff_y .- f_beat) .* 1e-6, SFM.nonneg(abs.(pp_y))), yscale=:log10, label=L"PSD of Signal $y(t_k)$")
		Plots.hline!(p2, [mn_noise_function(0)], label="Measurement Noise Floor", linewidth=2.5, linestyle=:dash)
		Plots.plot!(p2, xlabel="Frequency / MHz",
					ylabel="Signal PSD / (a.u.)",
				   # xticks=((-3:3),[L"f_{\mathrm{AOFS}} %$(i < 0 ? '-' : '+') %$(abs(i))" for i in (-3:3)])
					)
		p2
	end
end

# ╔═╡ 465e237d-a655-4f7d-85db-63edda1efacb
# DSP starts here


# ╔═╡ d9bb7605-acde-4fb9-8a46-d94980be5e88
begin
	if run_all
		f_bw = 2 * f_beat * 0.98
		phase_raw = SFM.bandpass_and_phase_estimation(y, f_sampling, f_beat, f_bw);
		phase, δ_freq = SFM.detrend_phase(phase_raw, f_sampling);
	end
end;

# ╔═╡ 487db32b-66d7-4951-8f7f-506c63f83481
begin
	f_sampling_raw = f_sampling;
	use_mP_for_train = training_estimator == "S_INV - P" ? true : false;

	l_delay_est = l_delay_m; # m
	tau_est = l_delay_est / (3e8/1.5);
	ONLY_FIRST_N_PEAKS = 5000;
	n_avg_param_est = 250;
	train_mask_filter_func = (f, p) -> (f .> 0.5/tau_est) .& (f .< f_bw/2) ;
	
	# this function additionally masks out spourious peaks
	# train_mask_filter_func = (f, p) -> (f .> 0.6e6) .& (f .< 3.0e7) .& (p .< (f.^(-2.0) .* 0.75e4))# for agilent laser (Parks ring)
	
	n_avg_krr_train = 150;
	n_avg = 50;
	T = Float64;
end;

# ╔═╡ 40368e5d-3b1c-4cb1-b683-2772973b4640
# estimate exp parameters

# ╔═╡ fab954c6-57fc-4057-a2cf-20acf188c053
begin
	if run_all	
		ff_phase, pp_phase = SFM.welch(phase, f_sampling, n_avg)
		init_exp_params = SFM.estimate_experimental_parameters(phase, f_sampling, f_bw, T(l_delay_est); n_psd_avg=n_avg_param_est, only_first_n_peaks=ONLY_FIRST_N_PEAKS);
		
		# optional: improve estimates using fringe fitting
		exp_params = SFM.improve_experimental_parameters_fringe_fit(phase, f_sampling, init_exp_params; n_psd_avg=100, median_filter_width=25);
	end
end;

# ╔═╡ f78e8896-beaa-437a-9dce-8a89bf364a48
if run_all
	if use_mP_for_train
		f_train, p_train, S_xx_inv_train, idx_train, mask_train, krr_f_range_ = SFM.get_train_samples_P(phase, f_sampling, exp_params, n_psd_avg=n_avg_krr_train, min_snr_dB=min_snr_dB, N_TRAIN_MAX=krr_train_samples);
	else
		f_train, p_train, S_xx_inv_train, idx_train, mask_train, krr_f_range_ = SFM.get_train_samples(phase, f_sampling, exp_params, n_psd_avg=n_avg_krr_train, min_snr_dB=min_snr_dB, N_TRAIN_MAX=krr_train_samples);
	end
end;

# ╔═╡ db246e44-7544-40f5-8dae-fbbdc81b5d52
begin	
	if run_all
		Plots.gr()
		p3 = Plots.plot(yscale=:log10, xlabel="Frequency / MHz", ylabel=L"PN-PSD / $rad^2\;Hz^{-1}$", palette = palette)
		Plots.plot!(p3, SFM.ds(Array(exp_params.f_psd)[2:end] .* 1e-6, exp_params.p_psd[2:end]),
				   label=L"S_{\Delta\phi}(f)")
	
		# # S_INV
		Plots.plot!(p3, SFM.ds(
		# 	S_freqs,
			ff_phase .* 1e-6,
			SFM.calc_S_inv(ff_phase, pp_phase, exp_params)
		), label=L"$S_{\phi, INV} (f)$")
	
		# S_INV - P
		Plots.plot!(p3, SFM.ds(
		# 	S_freqs,
			ff_phase .* 1e-6,
			SFM.calc_S_inv_mP(ff_phase, pp_phase, exp_params)
		), label=L"$S_{\phi, INV} (f) - P(f)$")
		
		Plots.scatter!(p3, exp_params.fringes.f .* 1e-6, exp_params.fringes.p, label="Detected Peaks")	
		
		Plots.hline!(p3, [exp_params.S_eta.val], label=L"Estimated Noise Floor $S_\epsilon$", linewidth=2, color=:black, linestyle=:dash)

		# Plots.hline!(p3, [2/(f_sampling * beat_note_snr_lin)], label="true Noise floor")
		
		Plots.plot!(p3, SFM.ds(Array(exp_params.f_psd .* 1e-6)[2:end], noise_func_pn(exp_params.f_psd[2:end])), label=L"True PN $S_{\phi, \mathrm{true}}(f)$", linewidth=2)
		
		Plots.xlims!(p3, 0, f_bw/2 * 1e-6)
		Plots.ylims!(p3, exp_params.S_eta.val/10, 1e-5)
		p3
	end
end

# ╔═╡ f5307a34-b25a-4f3a-9602-87d2f12f77cf


# ╔═╡ 9d4021ad-1b4d-4330-a687-f1bd3d6f9575


# ╔═╡ 5042b6d4-483b-4e17-ae3f-d1ee204f3ebd


# ╔═╡ 99c2143f-5790-4f6c-9682-ef4698990002
begin
	if run_all
		@warn "Using custom train sample filter"
		f_mask_train = train_mask_filter_func(f_train, p_train)
		f_train_masked = f_train[f_mask_train]
		p_train_masked = p_train[f_mask_train]
		
		
		# krr_f_range = krr_f_range_
		krr_f_range = [
			maximum([minimum(f_train_masked),krr_f_range_[1]]),
			minimum([maximum(f_train_masked), krr_f_range_[2]])
		]
		println("KRR range: ", krr_f_range)
		
		
		# mask away
		# f_mask = (f .> 0) .& (f .< krr_f_range[2])
		f_mask = (ff_phase .> 0) .& (ff_phase .< f_bw / 2)
		
		ff_phase_masked = ff_phase[f_mask]
		pp_phase_masked = pp_phase[f_mask]
		S_freqs = ff_phase_masked
		S_freqs_inv = ff_phase_masked
	end
end

# ╔═╡ 15f6b77b-6eea-45de-bb5d-44aed57fa132
begin
	if run_all
		using FLoops
		
		lambdas = T.(10.0 .^ range(-5.0, 0.0, length=n_krr_grid))
		if xlog
		    sigmas = T.(range(0.1, 10.0, length=n_krr_grid)) # for xlog = true
		else
		    sigmas = T.(10.0 .^ range(log10(5e4), log10(5e6), length=n_krr_grid))
		end
		
		results = zeros(T, length(lambdas), length(sigmas))
		results_std = zeros(T, length(lambdas), length(sigmas))
		
		x_val_groupsize = 25
		n_folds = 8
		
		for (i,lambda) in enumerate(lambdas)
		    @floop for (j,sigma) in enumerate(sigmas)
		        scores = SFM.grouped_cross_validation(
					f_train, p_train,
					SFM.krr_score_func $ (;lambda=lambda, kernel_func=kernel_func, kernel_params=(sigma=sigma,), xlog=xlog, ylog=ylog);
					folds=n_folds, groupsize=x_val_groupsize)
		        results[i,j] = mean(scores)
		        results_std[i,j] = std(scores)
		    end
		end
		
		# find optimal parameters via argmin
		min_idx = argmin(results)
		lambda_opt = lambdas[min_idx[1]]
		sigma_opt = sigmas[min_idx[2]]
		println("Optimal lambda: ", lambda_opt, " Optimal sigma: ", sigma_opt)
		
		
		krr_obj = SFM.krr_train(kernel_func, f_train, p_train; lambda=lambda_opt, kernel_params=(sigma=sigma_opt,), xlog=xlog, ylog=ylog)
		
		######################################
		
		mask_krr_range = (ff_phase_masked .> krr_f_range[1]) .& (ff_phase_masked .< krr_f_range[2])
		S_krr = SFM.krr_predict_batched(krr_obj, ff_phase_masked[mask_krr_range]; batch_size=10_000)
	
	
		# PSE
		S_joint = SFM.calc_S_inv(ff_phase_masked, pp_phase_masked, exp_params)
		# set S_xx_joint to KRR estimate where mask_krr_range is true
		S_joint[mask_krr_range] = S_krr
		S_joint[ff_phase_masked .> krr_f_range[2]] .= NaN
		snr = S_joint ./ exp_params.S_eta.val
		
		# define PSE filter
		H_w_2_inv_full = 1.0 ./ abs.(SFM.H_w_squared(ff_phase_masked, exp_params))
		G_PSE_2 = H_w_2_inv_full ./ (1.0 .+ H_w_2_inv_full .* (1.0 ./ snr));
		S_pse = pp_phase_masked .* G_PSE_2;
	end
end

# ╔═╡ 19b982ce-28ee-48c5-b45e-c962d3fc89f7
@suppress begin
	if run_all
		Plots.pythonplot()
		# show hyperparameter landscape
		Plots.heatmap(sigmas, lambdas, log10.(results), xlabel="Sigma", ylabel="Lambda", title="Hyperparameter landscape", color=:inferno, colorbar_title="log(Average MSE)", yscale=:log10, xscale=xlog ? :identity : :log10, vmin=-3, vmax=1)
		Plots.scatter!([sigma_opt], [lambda_opt], color=:lime, marker=:star4, lw=3, markersize=20, label="Optimal parameters")
		Plots.xlims!(sigmas[1], sigmas[end])
		Plots.ylims!(lambdas[1], lambdas[end])
	end
end

# ╔═╡ 42e826fc-dcb0-4057-9208-68009a818408
begin
	if run_all
		Plots.gr()
		p4 = Plots.plot(yscale=:log10, xlabel="Frequency / MHz", ylabel=L"PN-PSD / $rad^2\;Hz^{-1}$", palette = palette, leg=:topright)
		
		Plots.plot!(p4, SFM.ds(Array(exp_params.f_psd)[2:end] .* 1e-6, exp_params.p_psd[2:end]),
				   label=L"S_{\Delta\phi}(f)")
		Plots.plot!(p4, SFM.ds(
		# 	S_freqs,
			ff_phase .* 1e-6,
			SFM.calc_S_inv(ff_phase, pp_phase, exp_params)
		), label=L"$S_{\phi, INV} (f)$")
		
		Plots.plot!(p4, SFM.ds(ff_phase_masked.* 1e-6, S_pse)...,
					label=L"$S_{\phi, PSE} (f)$")
	
		# training samples
		n_sd_train_samples = 1
		Plots.scatter!(p4, f_train_masked[1:n_sd_train_samples:end] .* 1e-6, p_train_masked[1:n_sd_train_samples:end], label="Training Samples", alpha=0.4)
		
		# S KRR
		Plots.plot!(p4, SFM.ds(ff_phase_masked[mask_krr_range] .*1e-6, S_krr)...,
				   label=L"$S_{\phi, KRR} (f)$", linewidth=3)
	
		# S true
		Plots.plot!(p4, SFM.ds(Array(exp_params.f_psd .* 1e-6)[2:end], noise_func_pn(exp_params.f_psd[2:end])), label=L"True PN $S_{\phi, \mathrm{true}}(f)$", linewidth=2, linestyle=:dash)
	
	
	
		
			# segments
		# this is just to get the right freequencies, a bit hacky
		ff_train, _ = SFM.welch(phase, f_sampling, n_avg_krr_train)
		for seg in SFM.find_continuous_segments(.!BitVector(mask_train))
		    # println("Segment: ", ff_train[seg[1]]*1e-6, " - ", ff_train[seg[2]]*1e-6)
		    # f_lower = max(1e-5, ff_train[seg[1]]*1e-6)
		    # f_upper = min(10, ff_train[seg[2]]*1e-6)
		    f_lower = ff_train[seg[1]]*1e-6
		    f_upper = ff_train[seg[2]]*1e-6
		    # println("Segment: ", f_lower, " - ", f_upper)
		    Plots.vspan!(p4, [f_lower, f_upper],
						 color=:grey, alpha=0.3,
						 label="")
		end
		
		Plots.hline!(p4, [exp_params.S_eta.val], label=L"$S_\epsilon$", linewidth=2, color=:black, linestyle=:dash)
		Plots.xlims!(p4, 0, f_bw/2 * 1e-6)
		Plots.ylims!(p4, exp_params.S_eta.val/10, 1e-5)
	end
end

# ╔═╡ Cell order:
# ╟─0742953a-7f12-4a94-a42a-9433544bcfc1
# ╠═788fdd18-cf33-40f5-8845-9a35bed9c34d
# ╠═d09ffcda-1c5e-11f0-2a95-13c00027ecfa
# ╠═c27eafaf-cb46-4e90-9492-945221ebcca7
# ╠═89282ccd-b3c4-4587-af1f-96b6582b820a
# ╠═111f1904-5d2d-41ea-ad40-8a1b85dcf37b
# ╠═a6eee18b-7672-4154-88c8-fddb2dc5b6cf
# ╟─5107465f-5c73-41dd-af15-df87dc0b3d76
# ╟─df222550-e9ca-4bc0-bb01-67b6079c16ef
# ╟─a09fa9b8-8350-4951-a133-e50f7182db79
# ╟─5c394cd1-c50f-4048-ac08-c6d0e90b1851
# ╠═8da77ef1-ba12-4883-9f25-88fad120772a
# ╠═487db32b-66d7-4951-8f7f-506c63f83481
# ╠═369fa5fe-f2c1-4676-ab68-c3f6a8292656
# ╠═81ccb6d3-6608-49d1-a554-04a11a7839b0
# ╠═2b4666d8-16cb-46d0-bd48-89450f3e70ba
# ╟─646204c0-4bde-4f69-bfa9-3c481e03fb76
# ╠═465e237d-a655-4f7d-85db-63edda1efacb
# ╠═d9bb7605-acde-4fb9-8a46-d94980be5e88
# ╠═40368e5d-3b1c-4cb1-b683-2772973b4640
# ╠═fab954c6-57fc-4057-a2cf-20acf188c053
# ╠═f78e8896-beaa-437a-9dce-8a89bf364a48
# ╠═db246e44-7544-40f5-8dae-fbbdc81b5d52
# ╠═f5307a34-b25a-4f3a-9602-87d2f12f77cf
# ╠═9d4021ad-1b4d-4330-a687-f1bd3d6f9575
# ╠═5042b6d4-483b-4e17-ae3f-d1ee204f3ebd
# ╠═99c2143f-5790-4f6c-9682-ef4698990002
# ╠═15f6b77b-6eea-45de-bb5d-44aed57fa132
# ╠═19b982ce-28ee-48c5-b45e-c962d3fc89f7
# ╠═42e826fc-dcb0-4057-9208-68009a818408
