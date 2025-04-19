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
# ╠═╡ show_logs = false
begin
	# use the local environment and force the correct downloads
	import Pkg
	Pkg.activate(".")
	Pkg.resolve()
	# include("src/SDSHFringeMitigation.jl");
	using Suppressor
	import SDSHFringeMitigation; const SFM = SDSHFringeMitigation;
	
	# plotting
	import Plots
	Plots.default(dpi=300)
	Plots.pythonplot()

	# general use imports
	using LaTeXStrings, ColorSchemes, PlutoUI
	using PartialFunctions, Statistics, LinearAlgebra, DSP, FLoops
end;

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

# ╔═╡ 8dd41f4c-ad49-42ca-b5d0-d7e2b5e7f741
println("Installing required packages. This may take a few minutes...")

# ╔═╡ 5107465f-5c73-41dd-af15-df87dc0b3d76
md"""
# Load experimental data:
**Custom external cavity DFB laser**

Measured with a optical delay path of approx.  530 m

Choose the file: $(@bind whichfile Select(["ECDFB_10MHz.trc", "ECDFB_3MHz.trc"]))
"""

# ╔═╡ e2bed533-d289-44da-b9bb-f44cc2fde683
begin
	if whichfile == "ECDFB_10MHz.trc"
		f_beat = 10e6 # location of the beat note
		ONLY_FIRST_N_PEAKS = 20;
	elseif whichfile == "ECDFB_3MHz.trc"
		f_beat = 3e6 # location of the beat note
		ONLY_FIRST_N_PEAKS = 5;
	end
	
	f_bw = 0.96 * 2 * f_beat # filter bandwidth of bandpass filter around beat note
	l_delay_est = 530.0
	
	file = "exp_data/$(whichfile)"
	
	# load file
	t, y, f_sampling, T = SFM.load_data(file)
	
	# # # show the signal
	f_signal, p_signal = SFM.welch(y, f_sampling, 15)
	p1 = Plots.plot(SFM.ds(f_signal, p_signal),
			   yscale=:log10, xlabel="Frequency / Hz",
			   ylabel="PSD / (a.u.)", title="Measured Signal PSD",
			   legend=false)
	p1
end

# ╔═╡ 5c394cd1-c50f-4048-ac08-c6d0e90b1851
md"""
# KRR Parameters

Kernel Function: $(@bind kernel_func Select([SFM.rbf!, SFM.matern12!, SFM.matern32!, SFM.matern52!, SFM.matern72!, SFM.matern92!]))

SNR Threshold: $(@bind min_snr_dB NumberField(1.0:1.0:30.0, default=1)) dB

N training samples: $(@bind krr_train_samples NumberField(300:100:1000, default=600))

KRR x-axis log scaling: $(@bind xlog CheckBox(default=false))
KRR y-axis log scaling: $(@bind ylog CheckBox(default=true))

Estimator from which training samples are selected: $(@bind training_estimator Select(["S_INV - P", "S_INV"]))

KRR Hyperparameter Gridsearch size: $(@bind n_krr_grid NumberField(5:5:50, default=20))
"""

# ╔═╡ 8da77ef1-ba12-4883-9f25-88fad120772a
min_snr_dB; krr_train_samples; n_krr_grid; md"""
# Run everything?
#### Click the checkmark to compute everything!
#### ➡️ $(@bind run_all CheckBox(default=false)) ⬅️


The checkmark resets if phase noise parameters are changed, to save on computation.
"""

# ╔═╡ 487db32b-66d7-4951-8f7f-506c63f83481
begin
	n_si = 1.5
	c = 3e8

	use_mP_for_train = training_estimator == "S_INV - P" ? true : false;

	# l_delay_est = l_delay_m; # m
	tau_est = l_delay_est / (3e8/1.5);
	n_avg_param_est = 250;
	# train_mask_filter_func = (f, p) -> (f .> 0.5/tau_est) .& (f .< f_bw/2) ;
	train_mask_filter_func = (f, p) -> (f .> 0.5/tau_est) .& (f .< 8e6) ;
	
	# this function additionally masks out spourious peaks
	# train_mask_filter_func = (f, p) -> (f .> 0.6e6) .& (f .< 3.0e7) .& (p .< (f.^(-2.0) .* 0.75e4))# for agilent laser (Parks ring)
	
	n_avg_krr_train = 150;
	n_avg = 50;

	palette=:Set1_6
end;

# ╔═╡ d9bb7605-acde-4fb9-8a46-d94980be5e88
begin
	if run_all
		phase_raw = SFM.bandpass_and_phase_estimation(y, f_sampling, f_beat, f_bw);
		phase, δ_freq = SFM.detrend_phase(phase_raw, f_sampling);

		ff_phase, pp_phase = SFM.welch(phase, f_sampling, n_avg)
		init_exp_params = @suppress SFM.estimate_experimental_parameters(phase, f_sampling, f_bw, T(l_delay_est); n_psd_avg=n_avg_param_est, only_first_n_peaks=ONLY_FIRST_N_PEAKS);
		
		# optional: improve estimates using fringe fitting
		exp_params = SFM.improve_experimental_parameters_fringe_fit(phase, f_sampling, init_exp_params; n_psd_avg=100, median_filter_width=25);

		if use_mP_for_train
			f_train, p_train, S_xx_inv_train, idx_train, mask_train, krr_f_range_ = SFM.get_train_samples_P(phase, f_sampling, exp_params, n_psd_avg=n_avg_krr_train, min_snr_dB=min_snr_dB, N_TRAIN_MAX=krr_train_samples);
		else
			f_train, p_train, S_xx_inv_train, idx_train, mask_train, krr_f_range_ = SFM.get_train_samples(phase, f_sampling, exp_params, n_psd_avg=n_avg_krr_train, min_snr_dB=min_snr_dB, N_TRAIN_MAX=krr_train_samples);
		end
	end
end;

# ╔═╡ d70695fc-1c04-4473-8f48-176b4ee1a10a
md"""
## Conventional PN-PSD equalization methods
"""

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
		
		Plots.xlims!(p3, 0, f_bw/2 * 1e-6)
		Plots.ylims!(p3, exp_params.S_eta.val/10, 1e-5)
		p3
	end
end

# ╔═╡ 99c2143f-5790-4f6c-9682-ef4698990002
begin
	if run_all
		# @warn "Using custom train sample filter"
		f_mask_train = train_mask_filter_func(f_train, p_train)
		f_train_masked = f_train[f_mask_train]
		p_train_masked = p_train[f_mask_train]
		
		krr_f_range = [
			maximum([minimum(f_train_masked),krr_f_range_[1]]),
			minimum([maximum(f_train_masked), krr_f_range_[2]])
		]
		# println("KRR range: ", krr_f_range)
		
		# mask away
		# f_mask = (f .> 0) .& (f .< krr_f_range[2])
		f_mask = (ff_phase .> 0) .& (ff_phase .< f_bw / 2)
		
		ff_phase_masked = ff_phase[f_mask]
		pp_phase_masked = pp_phase[f_mask]
		S_freqs = ff_phase_masked
		S_freqs_inv = ff_phase_masked
	end
end;

# ╔═╡ 15f6b77b-6eea-45de-bb5d-44aed57fa132
begin
	if run_all
		# using FLoops
		
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
		# println("Optimal lambda: ", lambda_opt, " Optimal sigma: ", sigma_opt)
		
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
end;

# ╔═╡ 3dea087f-f796-4bf9-8d7d-4410bb007194
md"""
## KRR Hyperparameter grid search using grouped cross-validation
"""

# ╔═╡ 19b982ce-28ee-48c5-b45e-c962d3fc89f7
# ╠═╡ show_logs = false
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

# ╔═╡ b9564291-4cdb-456b-88db-4f70ec8f864f
md"""
## Data-driven PSE filter PN-PSE estimation
"""

# ╔═╡ 42e826fc-dcb0-4057-9208-68009a818408
begin
	if run_all
		Plots.gr()
		p5 = Plots.plot(yscale=:log10, xlabel="Frequency / MHz", ylabel=L"PN-PSD / $rad^2\;Hz^{-1}$", palette = palette, leg=:topright)
		
		Plots.plot!(p5, SFM.ds(Array(exp_params.f_psd)[2:end] .* 1e-6, exp_params.p_psd[2:end]),
				   label=L"S_{\Delta\phi}(f)")
		Plots.plot!(p5, SFM.ds(
		# 	S_freqs,
			ff_phase .* 1e-6,
			SFM.calc_S_inv(ff_phase, pp_phase, exp_params)
		), label=L"$S_{\phi, INV} (f)$")
		
		Plots.plot!(p5, SFM.ds(ff_phase_masked.* 1e-6, S_pse)...,
					label=L"$S_{\phi, PSE} (f)$")
	
		# training samples
		n_sd_train_samples = 5
		Plots.scatter!(p5, f_train_masked[1:n_sd_train_samples:end] .* 1e-6, p_train_masked[1:n_sd_train_samples:end], label="Training Samples", alpha=0.4)
		
		# S KRR
		Plots.plot!(p5, SFM.ds(ff_phase_masked[mask_krr_range] .*1e-6, S_krr)...,
				   label=L"$S_{\phi, KRR} (f)$", linewidth=3)
		
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
		    Plots.vspan!(p5, [f_lower, f_upper],
						 color=:grey, alpha=0.3,
						 label="")
		end
		
		Plots.hline!(p5, [exp_params.S_eta.val], label=L"$S_\epsilon$", linewidth=2, color=:black, linestyle=:dash)
		Plots.xlims!(p5, 0, f_bw/2 * 1e-6)
		Plots.ylims!(p5, exp_params.S_eta.val/10, 1e-5)
		p5
	end
end

# ╔═╡ 0d9e7ab1-fdee-4293-90b4-a8bff74f48be
begin
	if run_all
		Plots.gr()
		p4 = Plots.plot(xlabel="Frequency / Hz", ylabel=L"PN-PSD / $rad^2\;Hz^{-1}$", palette = palette, leg=:topright, xscale=:log10, yscale=:log10)
		
		Plots.plot!(p4, SFM.ds(Array(exp_params.f_psd)[2:end], exp_params.p_psd[2:end]),
				   label=L"S_{\Delta\phi}(f)")
		Plots.plot!(p4, SFM.ds(
		# 	S_freqs,
			ff_phase[2:end],
			SFM.calc_S_inv(ff_phase, pp_phase, exp_params)[2:end]
		), label=L"$S_{\phi, INV} (f)$")
		
		Plots.plot!(p4, SFM.ds(ff_phase_masked[2:end], S_pse[2:end])...,
					label=L"$S_{\phi, PSE} (f)$")
	
		# training samples
		Plots.scatter!(p4, f_train_masked[1:n_sd_train_samples:end], p_train_masked[1:n_sd_train_samples:end], label="Training Samples", alpha=0.4)
		
		# S KRR
		Plots.plot!(p4, SFM.ds(ff_phase_masked[mask_krr_range][2:end], S_krr[2:end])...,
				   label=L"$S_{\phi, KRR} (f)$", linewidth=3)
		
		for seg in SFM.find_continuous_segments(.!BitVector(mask_train))[2:end]
		    f_lower = ff_train[seg[1]]
		    f_upper = ff_train[seg[2]]
		    # println("Segment: ", f_lower, " - ", f_upper)
		    Plots.vspan!(p4, [f_lower, f_upper],
						 color=:grey, alpha=0.3,
						 label="")
		end
		
		Plots.hline!(p4, [exp_params.S_eta.val], label=L"$S_\epsilon$", linewidth=2, color=:black, linestyle=:dash)
		Plots.xlims!(p4, ff_phase[2], f_bw/2)
		Plots.ylims!(p4, exp_params.S_eta.val/10, maximum(S_pse))
		p4
	end
end

# ╔═╡ Cell order:
# ╟─0742953a-7f12-4a94-a42a-9433544bcfc1
# ╟─8dd41f4c-ad49-42ca-b5d0-d7e2b5e7f741
# ╟─d09ffcda-1c5e-11f0-2a95-13c00027ecfa
# ╟─5107465f-5c73-41dd-af15-df87dc0b3d76
# ╟─e2bed533-d289-44da-b9bb-f44cc2fde683
# ╟─5c394cd1-c50f-4048-ac08-c6d0e90b1851
# ╟─8da77ef1-ba12-4883-9f25-88fad120772a
# ╟─487db32b-66d7-4951-8f7f-506c63f83481
# ╟─d9bb7605-acde-4fb9-8a46-d94980be5e88
# ╟─d70695fc-1c04-4473-8f48-176b4ee1a10a
# ╟─db246e44-7544-40f5-8dae-fbbdc81b5d52
# ╟─99c2143f-5790-4f6c-9682-ef4698990002
# ╟─15f6b77b-6eea-45de-bb5d-44aed57fa132
# ╟─3dea087f-f796-4bf9-8d7d-4410bb007194
# ╟─19b982ce-28ee-48c5-b45e-c962d3fc89f7
# ╟─b9564291-4cdb-456b-88db-4f70ec8f864f
# ╟─42e826fc-dcb0-4057-9208-68009a818408
# ╟─0d9e7ab1-fdee-4293-90b4-a8bff74f48be
