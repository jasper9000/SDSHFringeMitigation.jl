using Measurements, Peaks, Serialization, LsqFit, Statistics
# import Plots

struct Fringes{T}
    f::AbstractArray{T, 1}
    p::AbstractArray{T, 1}
end

struct ExperimentalParameters{T}
    τ_d::Measurement{T}
    f_d::Measurement{T}
    l_d::Measurement{T}
    S_eta::Measurement{T}
    fringes::Fringes
    f_psd::AbstractArray{T, 1}
    p_psd::AbstractArray{T, 1}
end

# Define a method to save the struct to a file using Serialization
function save_exp_params(parameters::ExperimentalParameters, filename::String)
    open(filename, "w") do file
        serialize(file, parameters)
    end
end

# Define a method to load the struct from a file using Serialization
function load_exp_params(filename::String)::ExperimentalParameters
    open(filename, "r") do file
        return deserialize(file)
    end
end

function Base.show(io::IO, p::ExperimentalParameters{T}) where T
    println(io, "τ_d: \t", p.τ_d)
    println(io, "f_d: \t", p.f_d)
    println(io, "l_d: \t", p.l_d)
    println(io, "S_eta: \t", p.S_eta)
end

# override plot function for ExperimentalParameters
# function Plots.plot(p::ExperimentalParameters{T}; n_sig=2) where T
#     pp = Plots.plot(ds(Array(p.f_psd), nonneg(Array(p.p_psd)))...,label="Phase noise PSD", yscale=:log10, color=1)
#     Plots.scatter!(pp, p.fringes.f, p.fringes.p, label="Detected Peaks", color=2)
#     Plots.hline!(pp, [p.S_eta.val], label="", color=3, linewidth=2)
#     Plots.hspan!(pp, [max.(1e-14, p.S_eta.val - n_sig * p.S_eta.err), p.S_eta + n_sig * p.S_eta.err], label="Noise floor S_η ± $n_sig σ", alpha=0.15, color=3)
#     pp

#     # fringes
#     fringes_unc = (1:length(p.fringes.f)) * p.f_d
#     for f in fringes_unc
#         Plots.vline!([f.val], label="", color=4)
#         Plots.vspan!([f.val - n_sig * f.err, f.val + n_sig * f.err], label="", alpha=0.15, color=4)
#     end
#     # single label for fringes
#     Plots.vspan!([], label="Extracted fringes ± $n_sig σ", color=4)


#     # labels
#     Plots.plot!(xlabel="Frequency [Hz]", ylabel="Phase noise PSD [rad²/Hz]")
# end

function estimate_experimental_parameters(phase::ArrT, f_sampling::T, f_bw::T, l_delay_est::T; n_psd_avg::Int=100, only_first_n_peaks=nothing, s_eta_samples=4) where {T<:AbstractFloat, ArrT<:AbstractArray{T}}
    # calc heavily averaged PSD to estimate time delay and noise floor
    f, p = welch(phase, f_sampling, n_psd_avg; p_overlap=0.5, window=DSP.hanning)

    τ_d_estimated = l_delay_est / (3e8/1.5)

    # mask PSD to find peaks
    # f_mask_lower = min(1e5, 0.5/τ_d_estimated)
    f_mask_lower = 0.5/τ_d_estimated
    # println("Mask lower bound: ", f_mask_lower)
    f_mask = (f .> f_mask_lower) .& (f .< f_bw / 2)
    f_resolution = f[2] - f[1]
    # f_min_peak_dist = 1.0e5
    f_min_peak_dist = 0.25 / τ_d_estimated
    # println("Minimum peak distance: ", f_min_peak_dist)
    min_peak_dist = Int(floor(f_min_peak_dist / f_resolution))

    # on cpu
    peak_ids, peak_power = findminima(log10.(Array(p[f_mask])), min_peak_dist) |> peakproms!(;strict=false, min=0.2)
    f_peaks = f[f_mask][peak_ids]
    peak_power = 10.0.^(peak_power)
    # println(f"Found $(length(f_peaks)) peaks.")
    println("Found ", length(f_peaks), " peaks.")

    if !isnothing(only_first_n_peaks)
        println(only_first_n_peaks)
        # if it is an integer:
        if typeof(only_first_n_peaks) <: Int
            use_idx = 1:only_first_n_peaks
        # if iterable
        elseif typeof(only_first_n_peaks) <: AbstractArray
            use_idx = only_first_n_peaks
        else
            error("only_first_n_peaks must be an Int or an iterable.")
        end

        # check if use_idx is within bounds and simply use the peaks that are in the use_idx
        idx_max = findfirst(use_idx .> length(f_peaks))
        if !isnothing(idx_max)
            use_idx = use_idx[1:idx_max-1]
        end

        peak_ids = peak_ids[use_idx]
        f_peaks = f_peaks[use_idx]
        peak_power = peak_power[use_idx]
        println("Only using ", length(use_idx), " peaks.")
    end

    # f_0 = f_peaks[1]
    # f_diffs = abs.(f_peaks[2:end] .- f_0) ./ (1:length(f_peaks)-1)
    if length(f_peaks) > 1
        # f_diffs = abs.(diff(f_peaks))
        # f_diff = mean(f_diffs) ± std(f_diffs)

        # DIFFERENTT APPROADCH
        f_diffs = f_peaks ./ use_idx
        f_diff = mean(f_diffs) ± std(f_diffs)

    else
        f_diff = f_peaks[1] ± 0.0
    end

    t_diff = 1 / f_diff
    n_si = 1.5
    c = 3e8
    l_diff = t_diff .* T(c / n_si)


    # to estimate noise floor S_eta, we look at a few samples around the peak and take the median
    # s_eta_samples = 4
    s_eta_powers = Array{Float64}(undef, length(peak_ids))
    for (i, idx) in enumerate(peak_ids)
        s_eta_powers[i] = median(p[f_mask][max(1, idx-s_eta_samples÷2):min(end, idx+s_eta_samples÷2)])
    end
    # S_eta = median(s_eta_powers) ± std(s_eta_powers)

    # insetad, get the lower 10% quantile of the peak powers
    quant = quantile(s_eta_powers, 0.35)
    # get indices
    p10_powers = s_eta_powers[s_eta_powers .< quant]
    S_eta = median(p10_powers) ± std(p10_powers)

    # S_eta = mean(peak_power) ± std(peak_power)
    # use median instead of mean
    # S_eta = median(peak_power) ± std(peak_power)
    # make struct and make sure all arrays are of type Array
    f_peaks = Array(f_peaks)
    peak_power = Array(peak_power)
    params = ExperimentalParameters(
        t_diff,
        f_diff,
        l_diff,
        S_eta,
        # Fringes(f_peaks, peak_power),
        Fringes(f_peaks, s_eta_powers),
        Array(f[f_mask]),
        Array(p[f_mask])
    )
    return params
end

function improve_experimental_parameters_fringe_fit(phase::ArrT, f_sampling::T, exp_params::ExperimentalParameters; percent_fringe_fit_range=0.5::T, n_psd_avg::Int=4, median_filter_width::Int=100) where {T<:AbstractFloat, ArrT<:AbstractArray{T}}
    ff, psd = welch(phase, f_sampling, n_psd_avg; p_overlap=0.5, window=DSP.hanning)
    # fit model for a single fringe.
    # we assume here that the phase noise PSD is approx constant over one fringe
    model(f, p) = 10.0.^(p[1]) .* (2 .* (1 .- cos.(2pi * f .* p[2])) ) .+ 10.0.^p[3]
    # y log fit model for improved fitting
    model_ylog(f, p) = log10.(abs.(model(f, p)))
    
    f_diff_guess = Measurements.value(exp_params.f_d)
    f_log_factor = 10.0^floor(log10(f_diff_guess))
    # iterate over fringes
    fringe_fits = []
    for f_fringe in exp_params.fringes.f
        f_mask = (ff .> f_fringe - (percent_fringe_fit_range/2)*f_diff_guess) .& (ff .< f_fringe + (percent_fringe_fit_range/2)*f_diff_guess)
        
        xdata = ff[f_mask] / f_log_factor
        ydata = psd[f_mask]
        ydata_median = median_filter(ydata, median_filter_width)

        # use hanning window as fitting weights, such that the edges are weighted less
        ww = DSP.hanning(length(ff[f_mask])).^2

        p0 = [log10(ydata_median[1]), f_log_factor/f_diff_guess, log10(minimum(ydata_median))]
        p0_lower = [log10(1e-18), 0.5*f_log_factor/f_diff_guess, log10(1e-14)]
        p0_upper = [log10(1e-4), 1.5*f_log_factor/f_diff_guess, log10(1e-4)]

        # make sure the initial guess is within the bounds
        # println(p0)
        p0 = max.(p0, p0_lower)
        p0 = min.(p0, p0_upper)

        # linear model fit
        # fit = curve_fit(model, xdata, ydata_median, ww, p0, maxIter=1000,
        # lower=p0_lower, upper=p0_upper,
        # x_tol=1e-14, g_tol=1e-14)

        # y log model fit
        fit_ylog = curve_fit(model_ylog, xdata, log10.(abs.(ydata_median)), ww, p0; maxIter=1000,
        lower = p0_lower, upper = p0_upper, x_tol=1e-14, g_tol=1e-14
        )
        push!(fringe_fits, fit_ylog)
    end

    #### Use fits to update parameters ####
    taus = [fringe_fits[i].param[2] ±  mean(abs.(confidence_interval(fringe_fits[i])[2] .- fringe_fits[i].param[2])) for i in eachindex(fringe_fits)] / f_log_factor
    wm_tau = weightedmean(taus)
    wm_f_diff = 1 / wm_tau

    # note that this is an approximation as the real uncertainties are asymmentric
    S_etas = [10.0^fringe_fits[i].param[3] ± mean(abs.(10.0.^confidence_interval(fringe_fits[i])[3] .- 10.0^fringe_fits[i].param[3])) for i in eachindex(fringe_fits)]
    wm_S_eta = weightedmean(S_etas)

    n_si = 1.5
    c = 3e8
    wm_l_diff = wm_tau * T(c / n_si)

    which_idxs = round.(Int, Measurements.value.(exp_params.fringes.f / exp_params.f_d))

    new_exp_params = ExperimentalParameters(
        wm_tau,
        wm_f_diff,
        wm_l_diff,
        wm_S_eta,
        Fringes(which_idxs ./ taus, S_etas),
        # exp_params.f_psd,
        ff,
        # exp_params.p_psd
        psd
    )
    return new_exp_params
end


function median_filter(signal::Vector{T}, filter_length::Int) where T <: Real
    n = length(signal)
    half_length = filter_length ÷ 2
    filtered_signal = similar(signal)
    
    for i in 1:n
        # Determine the window boundaries
        start_idx = max(1, i - half_length)
        end_idx = min(n, i + half_length)
        
        # Extract the window
        window = signal[start_idx:end_idx]
        
        # Compute the median of the window
        filtered_signal[i] = median(window)
    end
    
    return filtered_signal
end