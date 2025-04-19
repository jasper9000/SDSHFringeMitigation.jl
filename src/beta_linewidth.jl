using Measurements

"Definition of beta line function"
beta_line = f -> 8*log(2)*f/(pi^2)

"Calculate the beta linewidth from the power spectral density (PSD) and the frequency array."
function beta_lw(f_psd::AbstractArray, p_fn_psd::ArrT) where ArrT<:AbstractArray
    beta_mask = p_fn_psd .> beta_line.(f_psd);
    lws = sqrt.(8*log(2)*fast_cumsum(reverse(p_fn_psd[beta_mask])) * (f_psd[2] - f_psd[1]) );
    f_onset = reverse(f_psd[beta_mask])
    return f_onset, lws
end

"Calculate the beta linewidth from the power spectral density (PSD) and the frequency array. In this version, the beta linewidth is calculated for the full range of frequencies starting at the first frequency above the beta line. This version struggles with spectra that feature spurious freqeuncies."
function beta_lw_full(f_psd::AbstractArray, p_fn_psd::AbstractArray)
    beta_mask_max = findlast(Array(p_fn_psd .> beta_line.(f_psd)));
    beta_mask = falses(length(f_psd))
    beta_mask[findfirst(f_psd .> 0):beta_mask_max] .= true
    lws = sqrt.(8*log(2)*fast_cumsum(reverse(p_fn_psd[beta_mask])) * (f_psd[2] - f_psd[1]) );
    f_onset = reverse(f_psd[beta_mask])
    return f_onset, lws
end

"Alias for cumsum to allow seamless use of fast_cumsum."
function fast_cumsum(x::AbstractVector{<:Number})
    cumsum(x)
end

"Implementation of cumsum for Measurement types for speed."
function fast_cumsum(x::AbstractVector{<:Measurement})
    val = x[1].val
    u2 = x[1].err ^ 2

    vals = zero(Measurements.value.(x))
    u2s = zero(Measurements.value.(x))

    vals[1] = val
    u2s[1] = u2

    for i in 2:length(x)
        val += x[i].val
        u2 += x[i].err^2
        vals[i] = val
        u2s[i] = u2
    end
    vals .± sqrt.(u2s)
end