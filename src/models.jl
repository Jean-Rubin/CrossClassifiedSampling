Base.@kwdef struct PopParam
    N_R::Int64
    N_C::Int64
end

Base.@kwdef struct BaseModel
    μ::Float64
    σ_R::Float64
    σ_C::Float64
    σ_RC::Float64
end

Base.@kwdef struct RatioPoisson
    base::BaseModel
    p::Float64
end

Base.@kwdef struct RatioIndep
    base::BaseModel
    p::Float64
end

Base.@kwdef struct RatioMixPoisson
    base::BaseModel
    α::Float64
    β::Float64
    μ_x::Float64
    μ_y::Float64
end

Base.@kwdef struct RatioMixSimple
    base::BaseModel
    base1::BaseModel
    base2::BaseModel
    α::Float64
    β::Float64
    μ_x::Float64
    μ_y::Float64
end


# --------------------------------------------------
# Model generation design
# --------------------------------------------------

function generate_base(model::BaseModel, param::PopParam)
    U = rand(Normal(0, model.σ_R), param.N_R, 1)
    V = rand(Normal(0, model.σ_C), 1, param.N_C)
    W = rand(Normal(0, model.σ_RC), param.N_R, param.N_C)

    Z = model.μ .+ U .+ V .+ W

    return Z
end

function generate(model::BaseModel, param::PopParam)
    Z = generate_base(model, param)

    Y = Array{Float64, 3}(undef, param.N_R, param.N_C, 1)

    @. Y[:, :, 1] = Z

    return Y
end

# Two variable model for ratio estimator
function generate(model::RatioPoisson, param::PopParam)
    Z = generate_base(model.base, param)

    Y = Array{Float64, 3}(undef, param.N_R, param.N_C, 2)

    @. Y[:, :, 1] = rand(Poisson(Z))
    @. Y[:, :, 2] = rand(Binomial(Y[:, :, 1], model.p))

    return Y
end

# Two variable model for ratio estimator
function generate(model::RatioIndep, param::PopParam)
    base₁ = model.base
    Z₁ = generate_base(base₁, param)
    base₂ = BaseModel(model.p * base₁.μ, base₁.σ_R, base₁.σ_C, base₁.σ_RC)
    Z₂ = generate_base(base₂, param)

    Y = Array{Float64, 3}(undef, param.N_R, param.N_C, 2)

    @. Y[:, :, 1] = Z₁
    @. Y[:, :, 2] = Z₂

    return Y
end

# Mixin two variable model for ratio estimator
function generate(model::RatioMixPoisson, param::PopParam)
    base = model.base
    Z₁ = generate_base(base, param)
    Z₂ = generate_base(base, param)
    Z₃ = generate_base(base, param)

    Y = Array{Float64, 3}(undef, param.N_R, param.N_C, 2)

    @. Y[:, :, 1] = rand(Poisson(model.μ_x + model.α * Z₁ + (1 - model.α) * Z₂))
    @. Y[:, :, 2] = rand(Poisson(model.μ_y + model.β * Z₁ + (1 - model.β) * Z₃))

    return Y
end

function generate(model::RatioMixSimple, param::PopParam)
    Z₁ = generate_base(model.base, param)
    Z₂ = generate_base(model.base1, param)
    Z₃ = generate_base(model.base2, param)

    Y = Array{Float64, 3}(undef, param.N_R, param.N_C, 2)

    @. Y[:, :, 1] = model.μ_x + model.α * Z₁ + (1 - model.α) * Z₂
    @. Y[:, :, 2] = model.μ_y + model.β * Z₁ + (1 - model.β) * Z₃

    return Y
end

function ratio(t)
    return t[2] / t[1]
end

function value(t)
    return t[1]
end
