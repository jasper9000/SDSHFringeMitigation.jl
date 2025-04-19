@doc """
Author: Jasper Riebesehl, 2025
"""

"Calculates modulus squared of interferometer transfer function H."
function H_w_squared(ff::ArrT, τ_d::T) where {T<:AbstractFloat, ArrT<:AbstractArray}
    return 2(1 .- cos.(2.0f0*π*ff*τ_d))
end

"Calculates modulus squared of interferometer transfer function H."
function H_w_squared(ff::AbstractArray, exp_params::ExperimentalParameters)
    return H_w_squared(ff, exp_params.τ_d.val)
end

"Calculates INV estimate of phase noise spectrum."
function calc_S_inv(phase::ArrT, f_sampling::T, τ_d::T; n_psd_avg=100, p_overlap=0.5, window=DSP.hanning) where {T<:AbstractFloat, ArrT<:AbstractArray{T}}
    f_psd, p_psd = welch(phase, f_sampling, n_psd_avg; p_overlap=p_overlap, window=window)
    return calc_S_inv(f_psd, p_psd, τ_d)
end

"Calculates INV estimate of phase noise spectrum."
function calc_S_inv(f_psd::ArrT, p_psd::ArrT, τ_d::T) where {T<:AbstractFloat, ArrT<:AbstractArray{T}}
    return p_psd ./ abs.(H_w_squared(f_psd, τ_d))
end

"Calculates INV estimate of phase noise spectrum."
function calc_S_inv(f_psd::ArrT, p_psd::ArrT, exp_params::ExperimentalParameters) where ArrT<:AbstractArray
    return calc_S_inv(f_psd, p_psd, exp_params.τ_d.val)
end 

"Calculates INV - P estimate of phase noise spectrum."
function calc_S_inv_mP(f_psd::ArrT, p_psd::ArrT, exp_params::ExperimentalParameters) where ArrT<:AbstractArray
    P = exp_params.S_eta.val ./ H_w_squared(f_psd, exp_params.τ_d.val)
    return calc_S_inv(f_psd, p_psd, exp_params.τ_d.val) .- P
end 

# "Calculates PSE estimate of phase noise spectrum."
# function calc_S_pse(f_psd::ArrT, p_psd::ArrT, krr_f_range::NTuple{2, T}, exp_params::ExperimentalParameters) where {ArrT<:AbstractArray, T<:Real}
#     S_joint = calc_S_inv(f_psd, p_psd, exp_params)
#     # set S_xx_joint to KRR estimate where mask_krr_range is true
#     mask_krr_range = (f .> krr_f_range[1]) .& (f .< krr_f_range[2])
#     S_joint[mask_krr_range] = S_krr
#     snr = S_joint ./ exp_params.S_eta.val

#     # define PSE filter
#     H_w_2_inv_full = 1.0 ./ abs.(H_w_squared(f_psd, exp_params))
#     G_PSE_2 = H_w_2_inv_full ./ (1.0 .+ H_w_2_inv_full .* (1.0 ./ snr));
#     # apply PSE filter
#     S_pse = p .* G_PSE_2;
#     return S_pse
# end

# "Calculates PSE estimate of phase noise spectrum."
# function calc_S_pse(phase::ArrT, f_sampling::T, krr_f_range::NTuple{2, T}, exp_params::ExperimentalParameters; n_psd_avg=100, p_overlap=0.5, window=DSP.hanning) where {ArrT<:AbstractArray, T<:Real}
#     f_psd, p_psd = welch(phase, f_sampling, n_psd_avg; p_overlap=p_overlap, window=window)
#     return calc_S_pse(f_psd, p_psd, krr_f_range, exp_params)
# end
