using CCS: BaseModel, PopParam, RatioMixSimple
using CCS: generate, simulate_completely, simulate_completely_var, value, ratio

function test_normal(σ_RC, n_R, n_C)
    @show σ_RC 
    @show n_R 
    @show n_C 
    Y = generate(
        BaseModel(μ = 200, σ_R = 5, σ_C = 5, σ_RC = σ_RC),
        PopParam(N_R = 1000, N_C = 1000)
    )

    nb_iter = 10000
    B = 1000

    simulate_completely_var(
    [Val(:Var1), Val(:Var2), Val(:Var3)],
        Y, value, n_R, n_C, nb_iter, B
    )

    simulate_completely(
        [Val(:Gross), Val(:RaoWu), Val(:Skinner)],
        [Val(:StdNormal), Val(:Percentile), Val(:RevPercentile)],
        Y, value, n_R, n_C, nb_iter, B
    )
end
test_normal(5, 5, 5)
test_normal(5, 10, 10)
test_normal(5, 10, 100)
test_normal(5, 100, 100)
test_normal(5, 500, 500)

test_normal(10, 5, 5)
test_normal(10, 10, 10)
test_normal(10, 10, 100)
test_normal(10, 100, 100)
test_normal(10, 500, 500)

test_normal(20, 5, 5)
test_normal(20, 10, 10)
test_normal(20, 10, 100)
test_normal(20, 100, 100)
test_normal(20, 500, 500)

function test_ratio_mix_simple(α, n_R, n_C)
    @show α 
    @show n_R 
    @show n_C 
    Y = generate(
        RatioMixSimple(
            base = BaseModel(μ = 0, σ_R = 5, σ_C = 5, σ_RC = 5),
            base1 = BaseModel(μ = 0, σ_R = 5, σ_C = 10, σ_RC = 15),
            base2 = BaseModel(μ = 0, σ_R = 5, σ_C = 10, σ_RC = 15),
            α = α,
            β = 0.5,
            μ_x = 100,
            μ_y = 300
        ),
        PopParam(N_R = 1000, N_C = 1000)
    )

    nb_iter = 10000
    B = 1000

    simulate_completely_var(
        [Val(:VarRatio1), Val(:VarRatio2), Val(:VarRatio3)],
        Y, ratio, n_R, n_C, nb_iter, B
    )

    simulate_completely(
        [Val(:Gross), Val(:RaoWu), Val(:Skinner)],
        [Val(:StdNormal), Val(:Percentile), Val(:RevPercentile)],
        Y, ratio, n_R, n_C, nb_iter, B
    )
end
test_ratio_mix_simple(0.1, 5, 5)
test_ratio_mix_simple(0.1, 10, 10)
test_ratio_mix_simple(0.1, 10, 100)
test_ratio_mix_simple(0.1, 100, 100)
test_ratio_mix_simple(0.1, 500, 500)

test_ratio_mix_simple(0.5, 5, 5)
test_ratio_mix_simple(0.5, 10, 10)
test_ratio_mix_simple(0.5, 10, 100)
test_ratio_mix_simple(0.5, 100, 100)
test_ratio_mix_simple(0.5, 500, 500)

test_ratio_mix_simple(0.8, 5, 5)
test_ratio_mix_simple(0.8, 10, 10)
test_ratio_mix_simple(0.8, 10, 100)
test_ratio_mix_simple(0.8, 100, 100)
test_ratio_mix_simple(0.8, 500, 500)

