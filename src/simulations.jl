function simulate_var(method::Val, Y, func::Function, n_R::Int64, n_C::Int64, nb_iter::Int64, B::Int64)
    N_R, N_C, _ = size(Y)
    boot_variance_results = Float64[]

    for iter in 1:nb_iter
        Y_sample = generate_sample(Y, n_R, n_C) 
        boot_variance = estimate_var(method, Y_sample, N_R, N_C, func, B)

        push!(boot_variance_results, boot_variance)
    end

    return boot_variance_results
end

function simulate_completely_var(methods_list, Y, func, n_R, n_C, nb_iter, B)
    print("Estimating true variance...\n")
    estimations = mcmc_true(Y, func, n_R, n_C, 1e5)
    true_variance = var(estimations)

    @show true_variance

    print("Bootstrapping...\n")
    for method in methods_list
        @show method
        boot_variance_results = simulate_var(method, Y, func, n_R, n_C, nb_iter, B)
        boot_variance_mean = mean(boot_variance_results)
        relative_bias = (boot_variance_mean - true_variance) / true_variance
        relative_stability = sqrt(mean((boot_variance_results .- true_variance).^2)) / true_variance
        @show relative_bias
        @show relative_stability
        print("\n")
    end

    return nothing
end


function simulate_bootstrap(method::Val, ci_methods, Y, func::Function, true_estimation::Float64, n_R::Int64, n_C::Int64, nb_iter::Int64, B::Int64)
    N_R, N_C, _ = size(Y)
    boot_variance_results = Float64[]
    coverage = fill(0., size(ci_methods))

    for iter in 1:nb_iter
        Y_sample = generate_sample(Y, n_R, n_C) 

        boot_estimations = generate_bootstrap(method, Y_sample, N_R, N_C, func, B)
        boot_variance = var(boot_estimations)

        for (i, ci_method) in enumerate(ci_methods)
            lower_ci, upper_ci = estimate_ci(ci_method, boot_estimations, 0.95)
            if lower_ci < true_estimation && true_estimation < upper_ci
                coverage[i] += 1.
            end
        end

        push!(boot_variance_results, boot_variance)
    end
    coverage ./= nb_iter

    return boot_variance_results, coverage
end


function simulate_completely(boot_methods, ci_methods, Y, func, n_R, n_C, nb_iter, B)
    print("Estimating true variance...\n")
    estimations = mcmc_true(Y, func, n_R, n_C, 1e5)
    true_variance = var(estimations)
    true_estimation = mean(estimations)

    @show true_variance

    print("Bootstrapping...\n")
    for method in boot_methods
        @show method
        boot_variance_results, coverage = simulate_bootstrap(method, ci_methods, Y, func, true_estimation, n_R, n_C, nb_iter, B)
        @show coverage
        boot_variance_mean = mean(boot_variance_results)
        relative_bias = (boot_variance_mean - true_variance) / true_variance
        relative_stability = sqrt(mean((boot_variance_results .- true_variance).^2)) / true_variance
        @show relative_bias
        @show relative_stability
        print("\n")
    end

    return nothing
end
