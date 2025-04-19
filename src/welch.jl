@doc """
Author: Jasper Riebesehl, 2025
"""

using DSP, FFTW

# make an enum for detrending options
@enum DetrendOption begin
    none
    linear
end

"CUDA- and CPU-compatible Welch method for estimating the power spectral density (PSD) of a signal."
function welch(signal::ArrT, f_sampling, n_avg; p_overlap=0.5, window=DSP.hanning, detrend::DetrendOption=linear) where ArrT<:AbstractArray
    N = length(signal)
    N_overlap = floor(Int, N * p_overlap) ÷ n_avg
    N_fft = floor(Int, N + (n_avg-1)*N_overlap) ÷ n_avg
    win = convert(ArrT, window(N_fft))
    segments = similar(signal, N_fft, n_avg)
    for i in 1:n_avg
        start = 1 + (i - 1) * (N_fft - N_overlap)
        seg = signal[start : start + N_fft-1]
        if detrend == linear
            seg, _ = linear_detrend(seg)
        end
        seg .*= win
        segments[:, i] .= seg
    end

    # apply rfft to the segments, plan it first
    F = plan_rfft(segments, 1)
    ft = F * segments

    # window function normalization
    ww = sum(win.^2)
    
    p = sum(2/(f_sampling*ww) * abs2.(ft), dims=2) / (n_avg)
    f = ArrT(rfftfreq(N_fft, f_sampling))
    return f, dropdims(p, dims=2)
end

"Simple function to perform a linear detrend."
function linear_detrend(data::ArrT; max_samples::Int = length(data)) where {ArrT <: AbstractArray}
    N = length(data)
    # Ensure max_samples does not exceed the length of the data
    sample_indices = if max_samples < N
        round.(Int, range(1, stop=N, length=max_samples))
    else
        1:N
    end
    
    # Downsample data and time index if necessary
    sampled_data = data[sample_indices]
    sampled_N = length(sampled_data)
    
    # Create the design matrix for linear regression
    X = similar(data, sampled_N, 2)
    X[:, 1] .= 1
    X[:, 2] .= sample_indices

    # Solve the least squares problem to find the slope and intercept
    coeffs = X \ sampled_data

    full_X = similar(data, N, 2)
    full_X[:, 1] .= 1
    full_X[:, 2] .= 1:N
    trend = full_X * coeffs
    
    # Subtract the trend from the original data
    # also return the coefficients
    return data .- trend, coeffs
end


