function estimate_ci(method::Val{:StdNormal}, boot_estimations, α)
    estimation = mean(boot_estimations)
    V_boot = var(boot_estimations)
    z = quantile(Normal(), (1 + α) / 2)

    lower = estimation - z * sqrt(V_boot)
    upper = estimation + z * sqrt(V_boot)

    return lower, upper
end


function estimate_ci(method::Val{:Percentile}, boot_estimations, α)
    frac = (1 - α) / 2

    sorted_estimations = sort(boot_estimations)
    B = length(boot_estimations)

    lower = sorted_estimations[floor(Int, frac * B)]
    upper = sorted_estimations[floor(Int, (1 - frac) * B)]

    return lower, upper
end


function estimate_ci(method::Val{:RevPercentile}, boot_estimations, α)
    frac = (1 - α) / 2

    estimation = mean(boot_estimations)
    sorted_estimations = sort(boot_estimations)
    B = length(boot_estimations)

    lower = 2 * estimation - sorted_estimations[floor(Int, (1 - frac) * B)]
    upper = 2 * estimation - sorted_estimations[floor(Int, frac * B)]

    return lower, upper
end
