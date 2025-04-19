using FFTW, Random

"Generates noise with custom power spectral density (PSD) using the inverse FFT method."
function generate_custom_noise(L::Int, f_sampling::Float64, noise_func::Function)
    rfreqs = rfftfreq(L+4, f_sampling)[2:end]
    n_rfreqs = length(rfreqs)
    empiric_psd_factor = sqrt(2 / (f_sampling * L))

    r_1 = (1/sqrt(2)) .* randn(n_rfreqs) .* sqrt.(noise_func(Array(rfreqs)))
    r_2 = (1/sqrt(2)) .* randn(n_rfreqs) .* sqrt.(noise_func(Array(rfreqs)))
    noise = r_1 .+ im .* r_2
    return real.(irfft(noise, L+2))[1:L] ./ empiric_psd_factor
end

"Custom noise function for a Lorentzian phase noise."
function noise_func_wiener(f, linewidth)
    return linewidth .* f.^(-2) ./ π
end

"Custom noise function for a white Gaussian noise."
function noise_func_white(f, noise_level_db)
    return 10.0 .^ (noise_level_db / 10.0)
end
