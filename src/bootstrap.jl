# --------------------------------------------------
# Bootstrap moments computation
# --------------------------------------------------

# Gross method
function generate_bootstrap(method::Val{:Gross}, Y_sample, N_R, N_C, func::Function, B::Int64)
    n_R, n_C, _ = size(Y_sample)
    boot_Y_pop = @view Y_sample[repeat(1:n_R, N_R ÷ n_R), repeat(1:n_C, N_C ÷ n_C), :]

    boot_estimations = Float64[]
    for boot_iter in 1:B
        boot_Y = generate_sample(boot_Y_pop, n_R, n_C)
        boot_estimation = estimate(boot_Y, func, N_R, N_C)

        push!(boot_estimations, boot_estimation)
    end

    return boot_estimations
end

function estimate_var(method::Val{:Gross}, Y_sample, N_R, N_C, func::Function, B::Int64)
    return var(generate_bootstrap(method, Y_sample, N_R, N_C, func, B))
end

# Rao and Wu method
function generate_bootstrap(method::Val{:RaoWu}, Y_sample, N_R, N_C, func::Function, B::Int64)
    n_R, n_C, _ = size(Y_sample)

    boot_estimations = Float64[]
    boot_Y = similar(Y_sample)
    m_R = Array{Int64}(undef, n_R, 1, 1)
    m_C = Array{Int64}(undef, 1, n_C, 1)

    for boot_iter in 1:B
        m_R[:, 1, 1] = rand(Multinomial(n_R - 1, repeat([1/n_R], n_R)))
        a_R = 1 .+ sqrt(1 - n_R / N_R) * (n_R * m_R / (n_R - 1) .- 1)
        m_C[1, :, 1] = rand(Multinomial(n_C - 1, repeat([1/n_C], n_C)))
        a_C = 1 .+ sqrt(1 - n_C / N_C) * (n_C * m_C / (n_C - 1) .- 1)

        @. boot_Y = a_R * a_C * Y_sample
        boot_estimation = estimate(boot_Y, func, N_R, N_C)

        push!(boot_estimations, boot_estimation)
    end

    return boot_estimations
end

function estimate_var(method::Val{:RaoWu}, Y_sample, N_R, N_C, func::Function, B::Int64)
    return var(generate_bootstrap(method, Y_sample, N_R, N_C, func, B))
end

# Skinner
function generate_bootstrap(method::Val{:Skinner}, Y_sample, N_R, N_C, func::Function, B::Int64)
    n_R, n_C, d = size(Y_sample)

    boot_estimations = Float64[]

    for boot_iter in 1:B
        boot_Y = @view Y_sample[rand(1:n_R, n_R - 1), rand(1:n_C, n_C - 1), :]
        boot_estimation = estimate(boot_Y, func, N_R, N_C)

        push!(boot_estimations, boot_estimation)
    end

    return boot_estimations
end

function estimate_var(method::Val{:Skinner}, Y_sample, N_R, N_C, func::Function, B::Int64)
    return var(generate_bootstrap(method, Y_sample, N_R, N_C, func, B))
end

# Variance estimator
# Note that func is assumed to be a projection
function estimate_var(method::Val{:Var1}, Y_sample, N_R, N_C, func::Function, B::Int64)
    n_R, n_C, d = size(Y_sample)
    f_R = n_R / N_R
    f_C = n_C / N_C

    Y_sample_mean_R = reshape(mean(Y_sample, dims = 2), n_R, 1, d)
    Y_sample_mean_C = reshape(mean(Y_sample, dims = 1), 1, n_C, d)
    Y_sample_mean_RC = reshape(mean(Y_sample, dims = (1,2)), 1, 1, d)

    V_hat_RC = sum((Y_sample .- Y_sample_mean_R .- Y_sample_mean_C .+ Y_sample_mean_RC).^2, dims = (1, 2)) / ((n_R - 1) * (n_C - 1))
    V_hat_R = sum((Y_sample_mean_R .- Y_sample_mean_RC).^2, dims = (1, 2)) / (n_R - 1) - (1 - f_C) * V_hat_RC / n_C
    V_hat_C = sum((Y_sample_mean_C .- Y_sample_mean_RC).^2, dims = (1, 2)) / (n_C - 1) - (1 - f_R) * V_hat_RC / n_R

    V_hat_all = (N_R * N_C)^2 * (
        (1 - f_R) * V_hat_R / n_R +
        (1 - f_C) * V_hat_C / n_C +
        (1 - f_R) * (1 - f_C) * V_hat_RC / (n_R * n_C)
    )

    return func(V_hat_all)
end

# Variance estimator but without V_RC for the final computation
# Note that func is assumed to be a projection
function estimate_var(method::Val{:Var2}, Y_sample, N_R, N_C, func::Function, B::Int64)
    n_R, n_C, d = size(Y_sample)
    f_R = n_R / N_R
    f_C = n_C / N_C

    Y_sample_mean_R = reshape(mean(Y_sample, dims = 2), n_R, 1, d)
    Y_sample_mean_C = reshape(mean(Y_sample, dims = 1), 1, n_C, d)
    Y_sample_mean_RC = reshape(mean(Y_sample, dims = (1,2)), 1, 1, d)

    V_hat_RC = sum((Y_sample .- Y_sample_mean_R .- Y_sample_mean_C .+ Y_sample_mean_RC).^2, dims = (1, 2)) / ((n_R - 1) * (n_C - 1))
    V_hat_R = sum((Y_sample_mean_R .- Y_sample_mean_RC).^2, dims = (1, 2)) / (n_R - 1) - (1 - f_C) * V_hat_RC / n_C
    V_hat_C = sum((Y_sample_mean_C .- Y_sample_mean_RC).^2, dims = (1, 2)) / (n_C - 1) - (1 - f_R) * V_hat_RC / n_R

    V_hat_all = (N_R * N_C)^2 * (
        (1 - f_R) * V_hat_R / n_R +
        (1 - f_C) * V_hat_C / n_C
    )

    return func(V_hat_all)
end

# Variance estimator without V_RC, and V_R and V_C are simplified
# Note that func is assumed to be a projection
function estimate_var(method::Val{:Var3}, Y_sample, N_R, N_C, func::Function, B::Int64)
    n_R, n_C, d = size(Y_sample)
    f_R = n_R / N_R
    f_C = n_C / N_C

    Y_sample_mean_R = reshape(mean(Y_sample, dims = 2), n_R, 1, d)
    Y_sample_mean_C = reshape(mean(Y_sample, dims = 1), 1, n_C, d)
    Y_sample_mean_RC = reshape(mean(Y_sample, dims = (1,2)), 1, 1, d)

    V_hat_R = sum((Y_sample_mean_R .- Y_sample_mean_RC).^2, dims = (1, 2)) / (n_R - 1)
    V_hat_C = sum((Y_sample_mean_C .- Y_sample_mean_RC).^2, dims = (1, 2)) / (n_C - 1)

    V_hat_all = (N_R * N_C)^2 * (
        (1 - f_R) * V_hat_R / n_R +
        (1 - f_C) * V_hat_C / n_C
    )

    return func(V_hat_all)
end

# Variance estimator for Ratio
# Thus assume that func is ratio
function estimate_var(method::Val{:VarRatio1}, Y_sample, N_R, N_C, func::Function, B::Int64)
    n_R, n_C, _ = size(Y_sample)
    t = estimate_t(Y_sample, N_R, N_C)
    Y_linear = similar(Y_sample, n_R, n_C, 1)

    Y_linear[:, :, 1] = (Y_sample[:, :, 2] - t[2]/t[1] * Y_sample[:, :, 1]) / t[1]

    return estimate_var(Val(:Var1), Y_linear, N_R, N_C, value, B)
end

# Assume that func is ratio
function estimate_var(method::Val{:VarRatio2}, Y_sample, N_R, N_C, func::Function, B::Int64)
    n_R, n_C, _ = size(Y_sample)
    t = estimate_t(Y_sample, N_R, N_C)
    Y_linear = similar(Y_sample, n_R, n_C, 1)

    Y_linear[:, :, 1] = (Y_sample[:, :, 2] - t[2]/t[1] * Y_sample[:, :, 1]) / t[1]

    return estimate_var(Val(:Var2), Y_linear, N_R, N_C, value, B)
end

# Assume that func is ratio
function estimate_var(method::Val{:VarRatio3}, Y_sample, N_R, N_C, func::Function, B::Int64)
    n_R, n_C, _ = size(Y_sample)
    t = estimate_t(Y_sample, N_R, N_C)
    Y_linear = similar(Y_sample, n_R, n_C, 1)

    Y_linear[:, :, 1] = (Y_sample[:, :, 2] - t[2]/t[1] * Y_sample[:, :, 1]) / t[1]

    return estimate_var(Val(:Var3), Y_linear, N_R, N_C, value, B)
end


# --------------------------------------------------
# 3 dimensional population
# --------------------------------------------------

# Gross method
function generate_bootstrap_3d(method::Val{:Gross}, Y_sample, N_1, N_2, N_3, func::Function, B::Int64)
    n_1, n_2, n_3, _ = size(Y_sample)
    boot_Y_pop = @view Y_sample[repeat(1:n_1, N_1 ÷ n_1), repeat(1:n_2, N_2 ÷ n_2), repeat(1:n_3, N_3 ÷ n_3), :]

    boot_estimations = Float64[]
    for boot_iter in 1:B
        boot_Y = generate_sample_3d(boot_Y_pop, n_1, n_2, n_3)
        boot_estimation = estimate_3d(boot_Y, func, N_1, N_2, N_3)

        push!(boot_estimations, boot_estimation)
    end

    return boot_estimations
end

function estimate_var_3d(method::Val{:Gross}, Y_sample, N_1, N_2, N_3, func::Function, B::Int64)
    return var(generate_bootstrap_3d(method, Y_sample, N_1, N_2, N_3, func, B))
end

# Rao and Wu method
function generate_bootstrap_3d(method::Val{:RaoWu}, Y_sample, N_1, N_2, N_3, func::Function, B::Int64)
    n_1, n_2, n_3, _ = size(Y_sample)

    boot_estimations = Float64[]
    boot_Y = similar(Y_sample)
    m_1 = Array{Int64}(undef, n_1, 1, 1, 1)
    m_2 = Array{Int64}(undef, 1, n_2, 1, 1)
    m_3 = Array{Int64}(undef, 1, 1, n_3, 1)

    for boot_iter in 1:B
        m_1[:, 1, 1, 1] = rand(Multinomial(n_1 - 1, repeat([1/n_1], n_1)))
        a_1 = 1 .+ sqrt(1 - n_1 / N_1) * (n_1 * m_1 / (n_1 - 1) .- 1)
        m_2[1, :, 1, 1] = rand(Multinomial(n_2 - 1, repeat([1/n_2], n_2)))
        a_2 = 1 .+ sqrt(1 - n_2 / N_2) * (n_2 * m_2 / (n_2 - 1) .- 1)
        m_3[1, 1, :, 1] = rand(Multinomial(n_3 - 1, repeat([1/n_3], n_3)))
        a_3 = 1 .+ sqrt(1 - n_3 / N_3) * (n_3 * m_3 / (n_3 - 1) .- 1)

        @. boot_Y = a_1 * a_2 * a_3 * Y_sample
        boot_estimation = estimate_3d(boot_Y, func, N_1, N_2, N_3)

        push!(boot_estimations, boot_estimation)
    end

    return boot_estimations
end

function estimate_var_3d(method::Val{:RaoWu}, Y_sample, N_1, N_2, N_3, func::Function, B::Int64)
    return var(generate_bootstrap_3d(method, Y_sample, N_1, N_2, N_3, func, B))
end

# Skinner
function generate_bootstrap_3d(method::Val{:Skinner}, Y_sample, N_1, N_2, N_3, func::Function, B::Int64)
    n_1, n_2, n_3, _ = size(Y_sample)

    boot_estimations = Float64[]

    for boot_iter in 1:B
        boot_Y = @view Y_sample[rand(1:n_1, n_1 - 1), rand(1:n_2, n_2 - 1), rand(1:n_3, n_3 - 1), :]
        boot_estimation = estimate_3d(boot_Y, func, N_1, N_2, N_3)

        push!(boot_estimations, boot_estimation)
    end

    return boot_estimations
end

function estimate_var_3d(method::Val{:Skinner}, Y_sample, N_1, N_2, N_3, func::Function, B::Int64)
    return var(generate_bootstrap_3d(method, Y_sample, N_1, N_2, N_3, func, B))
end

# Variance estimator
# Note that func is assumed to be a projection
function estimate_var_3d(method::Val{:Var1}, Y_sample, N_1, N_2, N_3, func::Function, B::Int64)
    n_1, n_2, n_3, _ = size(Y_sample)
    f_1 = n_1 / N_1
    f_2 = n_2 / N_2
    f_3 = n_3 / N_3

    Y_sample_mean_1 = mean(Y_sample, dims = (1))
    Y_sample_mean_2 = mean(Y_sample, dims = (2))
    Y_sample_mean_3 = mean(Y_sample, dims = (3))
    Y_sample_mean_12 = mean(Y_sample, dims = (1, 2))
    Y_sample_mean_13 = mean(Y_sample, dims = (1, 3))
    Y_sample_mean_23 = mean(Y_sample, dims = (2, 3))
    Y_sample_mean_123 = mean(Y_sample, dims = (1, 2, 3))

    V_hat_123 = sum((Y_sample .- Y_sample_mean_1 .- Y_sample_mean_2 .- Y_sample_mean_3 .+ Y_sample_mean_12 .+ Y_sample_mean_13 .+ Y_sample_mean_23 .- Y_sample_mean_123).^2, dims = (1, 2, 3)) / ((n_1 - 1) * (n_2 - 1) * (n_3 - 1))
    V_hat_12 = sum((Y_sample_mean_3 .- Y_sample_mean_13 .- Y_sample_mean_23 .+ Y_sample_mean_123).^2, dims = (1, 2, 3)) / ((n_1 - 1) * (n_2 - 1)) - (1 - f_3) * V_hat_123 / n_3
    V_hat_13 = sum((Y_sample_mean_2 .- Y_sample_mean_12 .- Y_sample_mean_23 .+ Y_sample_mean_123).^2, dims = (1, 2, 3)) / ((n_1 - 1) * (n_3 - 1)) - (1 - f_2) * V_hat_123 / n_2
    V_hat_23 = sum((Y_sample_mean_1 .- Y_sample_mean_12 .- Y_sample_mean_13 .+ Y_sample_mean_123).^2, dims = (1, 2, 3)) / ((n_2 - 1) * (n_3 - 1)) - (1 - f_1) * V_hat_123 / n_1
    V_hat_1 = sum((Y_sample_mean_23 .- Y_sample_mean_123).^2, dims = (1, 2, 3)) / (n_1 - 1) - (1 - f_2) * V_hat_12 / n_2 - (1 - f_3) * V_hat_13 / n_3 - (1 - f_2) * (1 - f_3) * V_hat_123 / (n_2 * n_3) 
    V_hat_2 = sum((Y_sample_mean_13 .- Y_sample_mean_123).^2, dims = (1, 2, 3)) / (n_2 - 1) - (1 - f_1) * V_hat_12 / n_1 - (1 - f_3) * V_hat_23 / n_3 - (1 - f_1) * (1 - f_3) * V_hat_123 / (n_1 * n_3)
    V_hat_3 = sum((Y_sample_mean_12 .- Y_sample_mean_123).^2, dims = (1, 2, 3)) / (n_3 - 1) - (1 - f_1) * V_hat_13 / n_1 - (1 - f_2) * V_hat_23 / n_2 - (1 - f_1) * (1 - f_2) * V_hat_123 / (n_1 * n_2)

    V_hat_all = (N_1 * N_2 * N_3)^2 * (
        (1 - f_1) * V_hat_1 / n_1 +
        (1 - f_2) * V_hat_2 / n_2 +
        (1 - f_3) * V_hat_3 / n_3 +
        (1 - f_1) * (1 - f_2) * V_hat_12 / (n_1 * n_2) +
        (1 - f_1) * (1 - f_3) * V_hat_13 / (n_1 * n_3) +
        (1 - f_2) * (1 - f_3) * V_hat_23 / (n_2 * n_3) +
        (1 - f_1) * (1 - f_2) * (1 - f_3) * V_hat_123 / (n_1 * n_2 * n_3)
    )

    return func(V_hat_all)
end

# Variance estimator but without V_12, V_13, V_23 and V_123 for the final computation
# Note that func is assumed to be a projection
function estimate_var_3d(method::Val{:Var2}, Y_sample, N_1, N_2, N_3, func::Function, B::Int64)
    n_1, n_2, n_3, _ = size(Y_sample)
    f_1 = n_1 / N_1
    f_2 = n_2 / N_2
    f_3 = n_3 / N_3

    Y_sample_mean_1 = mean(Y_sample, dims = (1))
    Y_sample_mean_2 = mean(Y_sample, dims = (2))
    Y_sample_mean_3 = mean(Y_sample, dims = (3))
    Y_sample_mean_12 = mean(Y_sample, dims = (1, 2))
    Y_sample_mean_13 = mean(Y_sample, dims = (1, 3))
    Y_sample_mean_23 = mean(Y_sample, dims = (2, 3))
    Y_sample_mean_123 = mean(Y_sample, dims = (1, 2, 3))

    V_hat_123 = sum((Y_sample .- Y_sample_mean_1 .- Y_sample_mean_2 .- Y_sample_mean_3 .+ Y_sample_mean_12 .+ Y_sample_mean_13 .+ Y_sample_mean_23 .- Y_sample_mean_123).^2, dims = (1, 2, 3)) / ((n_1 - 1) * (n_2 - 1) * (n_3 - 1))
    V_hat_12 = sum((Y_sample_mean_3 .- Y_sample_mean_13 .- Y_sample_mean_23 .+ Y_sample_mean_123).^2, dims = (1, 2, 3)) / ((n_1 - 1) * (n_2 - 1)) - (1 - f_3) * V_hat_123 / n_3
    V_hat_13 = sum((Y_sample_mean_2 .- Y_sample_mean_12 .- Y_sample_mean_23 .+ Y_sample_mean_123).^2, dims = (1, 2, 3)) / ((n_1 - 1) * (n_3 - 1)) - (1 - f_2) * V_hat_123 / n_2
    V_hat_23 = sum((Y_sample_mean_1 .- Y_sample_mean_12 .- Y_sample_mean_13 .+ Y_sample_mean_123).^2, dims = (1, 2, 3)) / ((n_2 - 1) * (n_3 - 1)) - (1 - f_1) * V_hat_123 / n_1
    V_hat_1 = sum((Y_sample_mean_23 .- Y_sample_mean_123).^2, dims = (1, 2, 3)) / (n_1 - 1) - (1 - f_2) * V_hat_12 / n_2 - (1 - f_3) * V_hat_13 / n_3 - (1 - f_2) * (1 - f_3) * V_hat_123 / (n_2 * n_3) 
    V_hat_2 = sum((Y_sample_mean_13 .- Y_sample_mean_123).^2, dims = (1, 2, 3)) / (n_2 - 1) - (1 - f_1) * V_hat_12 / n_1 - (1 - f_3) * V_hat_23 / n_3 - (1 - f_1) * (1 - f_3) * V_hat_123 / (n_1 * n_3)
    V_hat_3 = sum((Y_sample_mean_12 .- Y_sample_mean_123).^2, dims = (1, 2, 3)) / (n_3 - 1) - (1 - f_1) * V_hat_13 / n_1 - (1 - f_2) * V_hat_23 / n_2 - (1 - f_1) * (1 - f_2) * V_hat_123 / (n_1 * n_2)


    V_hat_all = (N_1 * N_2 * N_3)^2 * (
        (1 - f_1) * V_hat_1 / n_1 +
        (1 - f_2) * V_hat_2 / n_2 +
        (1 - f_3) * V_hat_3 / n_3
    )

    return func(V_hat_all)
end

# Variance estimator without V_12, V_13, V_23, and V_1, V_2, V_3 are simplified
# Note that func is assumed to be a projection
function estimate_var_3d(method::Val{:Var3}, Y_sample, N_1, N_2, N_3, func::Function, B::Int64)
    n_1, n_2, n_3, _ = size(Y_sample)
    f_1 = n_1 / N_1
    f_2 = n_2 / N_2
    f_3 = n_3 / N_3

    Y_sample_mean_12 = mean(Y_sample, dims = (1, 2))
    Y_sample_mean_13 = mean(Y_sample, dims = (1, 3))
    Y_sample_mean_23 = mean(Y_sample, dims = (2, 3))
    Y_sample_mean_123 = mean(Y_sample, dims = (1, 2, 3))

    V_hat_1_p = sum((Y_sample_mean_23 .- Y_sample_mean_123).^2, dims = (1, 2, 3)) / (n_1 - 1)
    V_hat_2_p = sum((Y_sample_mean_13 .- Y_sample_mean_123).^2, dims = (1, 2, 3)) / (n_2 - 1)
    V_hat_3_p = sum((Y_sample_mean_12 .- Y_sample_mean_123).^2, dims = (1, 2, 3)) / (n_3 - 1)

    V_hat_all = (N_1 * N_2 * N_3)^2 * (
        (1 - f_1) * V_hat_1_p / n_1 +
        (1 - f_2) * V_hat_2_p / n_2 +
        (1 - f_3) * V_hat_3_p / n_3
    )

    return func(V_hat_all)
end