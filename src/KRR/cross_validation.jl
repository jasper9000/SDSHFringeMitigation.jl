

function grouped_cross_validation(x_train::ArrT, y_train::ArrT, score_func::Function; folds=5, groupsize=50) where {T<:AbstractFloat, ArrT<:AbstractArray{T}}
    # created groups
    n_groups = length(x_train) ÷ groupsize + 1
    group_idxs = [[((ii-1)*groupsize+1):(ii*groupsize) for ii in 1:n_groups-1]; ((n_groups-1)*groupsize+1):length(x_train)]
    group_idxs = shuffleobs(group_idxs)

    # scores = []
    scores = zeros(T, folds)

    for (i, (idx_train, idx_val)) in enumerate(kfolds(collect(group_idxs), k=folds))
        train_idxs = sort(vcat(idx_train...))
        val_idxs = sort(vcat(idx_val...))
        x_train_fold = x_train[train_idxs]
        y_train_fold = y_train[train_idxs]
        x_val_fold = x_train[val_idxs]
        y_val_fold = y_train[val_idxs]
        
        score = score_func(x_train_fold, y_train_fold, x_val_fold, y_val_fold)
        # push!(scores, score)
        scores[i] = score
    end
    return scores
end

function krr_score_func(x_train, y_train, x_val, y_val; lambda=1e-3, kernel_func=rbf!, kernel_params=(sigma=1e5,), xlog=false, ylog=true)
    p_test = krr_fit_predict(kernel_func, x_train, y_train, x_val; lambda=lambda, kernel_params=kernel_params, xlog=xlog, ylog=ylog)
    if ylog # compare in log space to avoid crazy outliers
        return mean((log10.(p_test) - log10.(y_val)).^2)
    end
    return mean((p_test - y_val).^2)
end