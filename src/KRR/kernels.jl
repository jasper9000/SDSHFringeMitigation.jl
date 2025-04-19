using KernelAbstractions


@kernel function rbf!(K, x1, x2, @Const σ)
    i, j = @index(Global, NTuple)
    K[i,j] = exp(-(x1[i] - x2[j])^2 / (2σ^2))
end

@kernel function rbf_2!(K, x1, x2, @Const(σ), @Const(σ2))
    i, j = @index(Global, NTuple)
    K[i,j] = exp(-(x1[i] - x2[j])^2 / (2σ^2)) + exp(-(x1[i] - x2[j])^2 / (2σ2^2))
end

@kernel function matern12!(K, x1, x2, @Const σ)
    i, j = @index(Global, NTuple)
    K[i,j] = exp(-abs(x1[i] - x2[j]) / σ)
end

@kernel function matern32!(K, x1, x2, @Const σ)
    i, j = @index(Global, NTuple)
    d = abs(x1[i] - x2[j])
    sqrt_3 = 1.7320508f0
    K[i,j] = (1 + sqrt_3 * d / σ) * exp(-sqrt_3 * d / σ)
end

@kernel function matern52!(K, x1, x2, @Const σ)
    i, j = @index(Global, NTuple)
    d = abs(x1[i] - x2[j])
    sqrt_5 = 2.23606797749979f0
    K[i,j] = (1 + sqrt_5 * d / σ + 5/3 * (d / σ)^2) * exp(-sqrt_5 * d / σ)
end

@kernel function matern72!(K, x1, x2, @Const σ)
    i, j = @index(Global, NTuple)
    d = abs(x1[i] - x2[j])
    sqrt_7 = 2.64575131106459f0
    K[i,j] = (1 + sqrt_7 * d / σ + 7 * 2/5 * (d / σ)^2 + 7 * sqrt_7 * 1/15 * (d / σ)^3) * exp(-sqrt_7 * d / σ)
end

@kernel function matern92!(K, x1, x2, @Const σ)
    i, j = @index(Global, NTuple)
    d = abs(x1[i] - x2[j])
    sqrt_9 = 3.0f0
    K[i,j] = (
        1 + sqrt_9 * d / σ + 9 * 3/7* (d / σ)^2 + 9 * sqrt_9 * 2/21 * (d / σ)^3 + 9^2 * 1/105 * (d / σ)^4) * exp(-sqrt_9 * d / σ)
end
