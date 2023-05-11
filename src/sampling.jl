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
