# using KernelFunctions
using LinearAlgebra
using MLUtils

include("kernels.jl")

struct TrainedKRR{T<:AbstractFloat}
    x_train::AbstractArray{T}
    coefs::AbstractArray{T}
    lambda::T
    kernel_func::Function
    kernel_params::NamedTuple
    y_mean::T
    xlog::Bool
    ylog::Bool
end

function krr_train(kernel_func::Function, x_train::ArrT, y_train::ArrT; lambda=1e-3, kernel_params=nothing, xlog=false, ylog=true) where ArrT<:AbstractArray
    T = eltype(x_train)
    n_train = length(x_train)
    backend = get_backend(x_train)

    K = allocate(backend, T, n_train, n_train)

    if xlog
        x_train = log10.(x_train)
    end
    if kernel_func == rbf_2!
        kernel_func(get_backend(K))(K, x_train, x_train, kernel_params.sigma, kernel_params.sigma_2, ndrange=size(K))
    else
        kernel_func(get_backend(K))(K, x_train, x_train, kernel_params.sigma, ndrange=size(K))
    end 
    
    if ylog
        # coefs = (K + lambda * I) \ log10.(y_train)
        coefs = inv(K + lambda * I) * log10.(y_train)
        y_mean = log10(mean(y_train))
    else
        # coefs = (K + lambda * I) \ y_train
        coefs = inv(K + lambda * I) *  y_train
        y_mean = mean(y_train)
    end
    return TrainedKRR(x_train, coefs, lambda, kernel_func, kernel_params, y_mean, xlog, ylog)
end

function krr_predict(krr_obj::TrainedKRR, x_test::ArrT) where ArrT<:AbstractArray
    T = eltype(krr_obj.x_train)
    n_train = length(krr_obj.x_train)
    backend = get_backend(krr_obj.x_train)

    if krr_obj.xlog
        x_test = log10.(x_test)
    end

    Kstar = allocate(backend, T, length(x_test), n_train)
    krr_obj.kernel_func(backend)(Kstar, x_test, krr_obj.x_train, krr_obj.kernel_params.sigma, ndrange=size(Kstar))
    if krr_obj.ylog
        return 10.0 .^ (Kstar * krr_obj.coefs)
    else
        return Kstar * krr_obj.coefs #.+ krr_obj.y_mean
    end
end

function krr_predict_batched(krr_obj::TrainedKRR, x_test::AbstractArray{T}; batch_size::Int=10_000) where T<:AbstractFloat
    l = length(x_test)
    y_pred = allocate(get_backend(x_test), T, l)
    for i in 1:batch_size:l
        y_pred[i:min(i+batch_size, l)] .= krr_predict(krr_obj, x_test[i:min(i+batch_size, l)])
    end
    return y_pred #.+ krr_obj.y_mean
    # return 10.0 .^ (y_pred .+ krr_obj.y_mean)
end

function krr_fit_predict(kernel_func::Function, x_train::ArrT, y_train::ArrT, x_test::ArrT; lambda::T=1e-3, kernel_params=nothing, xlog=false, ylog=true) where {T<:AbstractFloat, ArrT<:AbstractArray{T}}
    k = krr_train(kernel_func, x_train, y_train; lambda=lambda, kernel_params=kernel_params, xlog=xlog, ylog=ylog)
    return krr_predict(k, x_test)
end

# function grouped_cross_validation(x_train::ArrT, y_train::ArrT, score_func::Function; folds=5, groupsize=50) where {T<:AbstractFloat, ArrT<:AbstractArray{T}}
#     # created groups
#     n_groups = length(x_train) ÷ groupsize + 1
#     group_idxs = [[((ii-1)*groupsize+1):(ii*groupsize) for ii in 1:n_groups-1]; ((n_groups-1)*groupsize+1):length(x_train)]
#     group_idxs = shuffleobs(group_idxs)

#     # scores = []
#     scores = zeros(T, folds)

#     for (i, (idx_train, idx_val)) in enumerate(kfolds(collect(group_idxs), k=folds))
#         train_idxs = sort(vcat(idx_train...))
#         val_idxs = sort(vcat(idx_val...))
#         x_train_fold = x_train[train_idxs]
#         y_train_fold = y_train[train_idxs]
#         x_val_fold = x_train[val_idxs]
#         y_val_fold = y_train[val_idxs]
        
#         score = score_func(x_train_fold, y_train_fold, x_val_fold, y_val_fold)
#         # push!(scores, score)
#         scores[i] = score
#     end
#     return scores
# end

# function krr_score_func(x_train, y_train, x_val, y_val; lambda=1e-3, kernel_func=rbf!, kernel_params=(sigma=1e5,), xlog=false, ylog=true)
#     p_test = krr_fit_predict(kernel_func, x_train, y_train, x_val; lambda=lambda, kernel_params=kernel_params, xlog=xlog, ylog=ylog)
#     if ylog # compare in log space to avoid crazy outliers
#         return mean((log10.(p_test) - log10.(y_val)).^2)
#     end
#     return mean((p_test - y_val).^2)
# end