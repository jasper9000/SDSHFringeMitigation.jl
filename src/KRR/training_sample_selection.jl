@doc """
Author: Jasper Riebesehl, 2025
"""

function get_train_samples(phase::ArrT, f_sampling::T, exp_params::ExperimentalParameters; n_psd_avg::Int=150, min_snr_dB::T=15.0, N_TRAIN_MAX::Int=2000, kwargs...) where {T<:AbstractFloat, ArrT<:AbstractArray}
    f_psd, p_psd = welch(phase, f_sampling, n_psd_avg; kwargs...)
    idx_train, mask_train, krr_f_range = get_train_idx(f_psd, p_psd, exp_params, min_snr_dB, N_TRAIN_MAX)
    S_xx_inv_train = calc_S_inv(f_psd, p_psd, exp_params.τ_d.val)
    f_train = f_psd[mask_train][idx_train]
    p_train = S_xx_inv_train[mask_train][idx_train]
    return f_train, p_train, S_xx_inv_train, idx_train, mask_train, krr_f_range
end

function get_train_samples(f_psd::ArrT, p_psd::ArrT, exp_params::ExperimentalParameters; min_snr_dB::T=15.0, N_TRAIN_MAX::Int=2000) where {T<:AbstractFloat, ArrT<:AbstractArray}
    idx_train, mask_train, krr_f_range = get_train_idx(f_psd, p_psd, exp_params, min_snr_dB, N_TRAIN_MAX)
    S_xx_inv_train = calc_S_inv(f_psd, p_psd, exp_params.τ_d.val)
    f_train = f_psd[mask_train][idx_train]
    p_train = S_xx_inv_train[mask_train][idx_train]
    return f_train, p_train, S_xx_inv_train, idx_train, mask_train, krr_f_range
end

function get_snr_mask_function(exp_params::ExperimentalParameters; min_snr_dB::T=15.0) where T<:AbstractFloat
    return (f_psd, p_psd) -> (10log10.(p_psd ./ exp_params.S_eta.val) .> min_snr_dB) .& (f_psd .> 1e3)
end

function get_train_idx(f_psd::ArrT, p_psd::ArrT, exp_params::ExperimentalParameters, min_snr_dB::T=15.0, N_TRAIN_MAX::Int=2000) where {T<:AbstractFloat, ArrT<:AbstractArray}
    mask_train = get_snr_mask_function(exp_params; min_snr_dB)(f_psd, p_psd)
    krr_f_range = [f_psd[findfirst(mask_train)], f_psd[findlast(mask_train)]]
    mask_train .&= (f_psd .> krr_f_range[1]) .& (f_psd .< krr_f_range[2])

    # select from training data
    idx_train = 1:(1 + sum(mask_train) ÷ N_TRAIN_MAX):sum(mask_train)
    f_train = f_psd[mask_train][idx_train]

    # update krr range to the selected range
    krr_f_range = [minimum(f_train), maximum(f_train)]
    return idx_train, mask_train, krr_f_range
end


function get_train_samples_P(phase::ArrT, f_sampling::T, exp_params::ExperimentalParameters; n_psd_avg::Int=150, min_snr_dB::T=15.0, N_TRAIN_MAX::Int=2000, kwargs...) where {T<:AbstractFloat, ArrT<:AbstractArray}
    f_psd, p_psd = welch(phase, f_sampling, n_psd_avg; kwargs...)
    idx_train, mask_train, krr_f_range = get_train_idx(f_psd, p_psd, exp_params, min_snr_dB, N_TRAIN_MAX)
    S_xx_inv_train = calc_S_inv(f_psd, p_psd, exp_params.τ_d.val)
    P_train = S_xx_inv_train .- exp_params.S_eta.val ./ H_w_squared(f_psd, exp_params.τ_d.val)
    f_train = f_psd[mask_train][idx_train]
    p_train = P_train[mask_train][idx_train]
    return f_train, p_train, P_train, idx_train, mask_train, krr_f_range
end