using CCS: BaseModel, PopParam, RatioMixSimple
using CCS: generate, simulate_completely, simulate_completely_var, value, ratio
using CCS: Model3D, PopParam3D
using CCS: simulate_completely_3d, simulate_completely_var_3d

function test_normal(σ_RC, n_R, n_C)
    @show σ_RC 
    @show n_R, n_C 
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

function test_normal_3d(σ_1, σ_12, σ_123, n_1, n_2, n_3)
    @show σ_1, σ_12, σ_123
    @show n_1, n_2, n_3 
    Y = generate(
        Model3D(μ = 200, σ_1 = σ_1, σ_2 = 5, σ_3 = 5, σ_12 = σ_12, σ_13 = 5, σ_23 = 5, σ_123 = σ_123),
        PopParam3D(N_1 = 1000, N_2 = 100, N_3 = 100)
    )

    nb_iter = 10000
    B = 1000

    simulate_completely_var_3d(
    [Val(:Var1), Val(:Var2), Val(:Var3)],
        Y, value, n_1, n_2, n_3, nb_iter, B
    )

    simulate_completely_3d(
        [Val(:Gross), Val(:RaoWu), Val(:Skinner)],
        [Val(:StdNormal), Val(:Percentile), Val(:RevPercentile)],
        Y, value, n_1, n_2, n_3, nb_iter, B
    )
end

open("res_3d.txt", "w") do out
    redirect_stdout(out) do
    @time test_normal_3d(5, 5, 5, 10, 10, 10)
    @time test_normal_3d(5, 5, 5, 10, 50, 50)
    @time test_normal_3d(5, 5, 5, 100, 10, 10)
    @time test_normal_3d(5, 5, 5, 100, 50, 50)
    @time test_normal_3d(5, 5, 5, 500, 10, 10)
    @time test_normal_3d(5, 5, 5, 500, 50, 50)
    @time test_normal_3d(5, 10, 5, 10, 10, 10)
    @time test_normal_3d(5, 10, 5, 10, 50, 50)
    @time test_normal_3d(5, 10, 5, 500, 10, 10)
    @time test_normal_3d(5, 10, 5, 500, 50, 50)
    @time test_normal_3d(5, 10, 5, 100, 10, 10)
    @time test_normal_3d(5, 10, 5, 100, 50, 50)
    @time test_normal_3d(5, 20, 5, 100, 10, 10)
    @time test_normal_3d(5, 20, 5, 100, 50, 50)
    @time test_normal_3d(5, 20, 5, 10, 10, 10)
    @time test_normal_3d(5, 20, 5, 10, 50, 50)
    @time test_normal_3d(5, 20, 5, 500, 10, 10)
    @time test_normal_3d(5, 20, 5, 500, 50, 50)
    @time test_normal_3d(5, 5, 10, 100, 10, 10)
    @time test_normal_3d(5, 5, 10, 100, 50, 50)
    @time test_normal_3d(5, 5, 10, 10, 10, 10)
    @time test_normal_3d(5, 5, 10, 10, 50, 50)
    @time test_normal_3d(5, 5, 10, 500, 10, 10)
    @time test_normal_3d(5, 5, 10, 500, 50, 50)
    @time test_normal_3d(5, 5, 20, 100, 10, 10)
    @time test_normal_3d(5, 5, 20, 100, 50, 50)
    @time test_normal_3d(5, 5, 20, 10, 10, 10)
    @time test_normal_3d(5, 5, 20, 10, 50, 50)
    @time test_normal_3d(5, 5, 20, 500, 10, 10)
    @time test_normal_3d(5, 5, 20, 500, 50, 50)
    end
end