# --------------------------------------------------
# Sampling and estimation helper
# --------------------------------------------------

function generate_sample(Y, n_R, n_C)
    N_R, N_C, _ = size(Y)
    S_R = sample(1:N_R, n_R, replace = false)
    S_C = sample(1:N_C, n_C, replace = false)

    return @view Y[S_R, S_C, :]
end

function estimate_t(Y_sample, N_R, N_C)
    n_R, n_C, _ = size(Y_sample) 
    
    return N_R * N_C / (n_R * n_C) * sum(Y_sample, dims = (1, 2))
end

function estimate(Y_sample, func, N_R, N_C)
    return func(estimate_t(Y_sample, N_R, N_C))
end


# 3 dimensional population
function generate_sample_3d(Y, n_1, n_2, n_3)
    N_1, N_2, N_3, _ = size(Y)
    S_1 = sample(1:N_1, n_1, replace = false)
    S_2 = sample(1:N_2, n_2, replace = false)
    S_3 = sample(1:N_3, n_3, replace = false)

    return @view Y[S_1, S_2, S_3, :]
end

function estimate_t_3d(Y_sample, N_1, N_2, N_3)
    n_1, n_2, n_3, _ = size(Y_sample)

    return N_1 * N_2 * N_3 / (n_1 * n_2 * n_3) * sum(Y_sample, dims = (1, 2, 3))
end

function estimate_3d(Y_sample, func, N_1, N_2, N_3)
    return func(estimate_t_3d(Y_sample, N_1, N_2, N_3))
end

# --------------------------------------------------
# True variance estimation from MC on estimations
# --------------------------------------------------

function mcmc_true(Y, func, n_R, n_C, nb_iter)
    N_R, N_C, _ = size(Y)
    estimations = Float64[]

    for iter in 1:nb_iter
        Y_sample = generate_sample(Y, n_R, n_C)
        estimation = estimate(Y_sample, func, N_R, N_C)

        push!(estimations, estimation)
    end

    return estimations
end

# 3 dimensional population
function mcmc_true_3d(Y, func, n_1, n_2, n_3, nb_iter)
    N_1, N_2, N_3, _ = size(Y)
    estimations = Float64[]

    for iter in 1:nb_iter
        Y_sample = generate_sample_3d(Y, n_1, n_2, n_3)
        estimation = estimate_3d(Y_sample, func, N_1, N_2, N_3)

        push!(estimations, estimation)
    end

    return estimations
end