using Adapt

function aliased_freq(f, f_sampling)
    return abs(f - round(f / f_sampling) * f_sampling)
end

function get_filter_band(f_center, f_bw, f_sampling)
    f_center_alias = aliased_freq(f_center, f_sampling)
    f_band = [maximum([1, f_center_alias - f_bw/2]), minimum([f_sampling / 2, f_center_alias + f_bw/2])]
    return f_band
end

function phase_estimation(signal::ArrT, f_sampling::T, f_beat::T) where {T<:AbstractFloat, ArrT<:AbstractArray{T}}
    ArrWrapper = typeof(signal).name.wrapper
    N = length(signal)

    # fourier transform
    ft = 2 * fft(signal)
    ft_freqs = fftfreq(N, f_sampling)

    # discrete Hilbert transform
    mask = ft_freqs .< 0
    ft[mask] .= 0

    # Inv Fourier transform, remove linear phase trend
    quadrature = ifft(ft) .* exp.(-im * 2 * pi * f_beat * ArrWrapper{T}(collect(0.0:N-1)) / f_sampling)

    # extract und unwrap phase
    b = angle.(quadrature)
    # do this on the CPU for now
    phase = unwrap(Array(b))

    # convert back to original array type
    return adapt(ArrWrapper, phase)
end

function phase_and_amp_estimation(signal::ArrT, f_sampling::T, f_beat::T) where {T<:AbstractFloat, ArrT<:AbstractArray{T}}
    ArrWrapper = typeof(signal).name.wrapper
    N = length(signal)

    # Fourier transform
    ft = 2 * fft(signal)
    ft_freqs = fftfreq(N, f_sampling)

    # discrete Hilbert transform
    mask = ft_freqs .< 0
    ft[mask] .= 0

    # Inv Fourier transform, remove linear phase trend
    quadrature = ifft(ft) .* exp.(-im * 2 * pi * f_beat * ArrWrapper{T}(collect(0.0:N-1)) / f_sampling)

    # extract und unwrap phase
    b = angle.(quadrature)
    # do this on the CPU for now
    phase = unwrap(Array(b))

    # extract amplitude (if needed)
    amp = abs.(quadrature)

    # convert back to original array type
    return adapt(ArrWrapper, phase), adapt(ArrWrapper, amp)
end


function bandpass_and_phase_estimation(signal::ArrT, f_sampling::T, f_beat::T, f_bw_filter::T) where {T<:AbstractFloat, ArrT<:AbstractArray{T}}
    ArrWrapper = typeof(signal).name.wrapper
    N = length(signal)

    f_band = get_filter_band(f_beat, f_bw_filter, f_sampling)

    # fourier transform
    ft = 2 * fft(signal)
    ft_freqs = fftfreq(N, f_sampling)

    # println("filter band: ", f_band)
    mask = (ft_freqs .>= f_band[1]) .& (ft_freqs .<= f_band[2])

    ft[.!mask] .= 0

    quadrature = ifft(ft) .* exp.(-im * 2 * pi * f_beat * ArrWrapper{T}(collect(0.0:N-1)) / f_sampling)

    b = angle.(quadrature)

    # do this on the CPU for now
    phase = unwrap(Array(b))

    # convert back to original array type
    return change_float_precision(adapt(ArrWrapper, phase), T)
end

"Wrapper function around linear_detrend."
function detrend_phase(phase::ArrT, f_sampling::T; n_linregress=1000) where {T<:AbstractFloat, ArrT<:AbstractArray{T}}
    ArrWrapper = typeof(phase).name.wrapper
    detrended_phase, coefs = linear_detrend(phase; max_samples=n_linregress)
    delta_f = f_sampling * coefs[2] / (2 * pi)

    return adapt(ArrWrapper, T.(detrended_phase)), delta_f
end

# function detrend_phase(phase::ArrT, f_sampling::T; n_linregress=1000) where {T<:AbstractFloat, ArrT<:AbstractArray{T}}
#     ArrWrapper = typeof(phase).name.wrapper
#     N = length(phase)
#     t = adapt(ArrWrapper, collect(0.0:N-1) / f_sampling)

#     # eqally spaced index of elements used for linear regression
#     idx = 1:N÷n_linregress:N

#     # do this on CPU
#     model = GLM.lm(@formula(phase ~ t),
#         # DataFrame(phase=convert.(Float64, phase[idx]), t=convert.(Float64, t[idx]))
#         DataFrame(phase=Array{Float64}(phase[idx]), t=Array{Float64}(t[idx]))
#         )

#     # return detrended phase and frequency
#     return adapt(ArrWrapper, T.(phase - GLM.coef(model)[2] * t)), GLM.coef(model)[2] / (2 * pi)
# end

