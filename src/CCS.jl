module CCS

using Distributions
using StatsBase

export PopParam, BaseModel, RatioPoisson, RatioIndep, RatioMixPoisson, RatioMixSimple
export generate, simulate_completely, simulate_completely_var

include("models.jl")
include("sampling.jl")
include("bootstrap.jl")
include("intervals.jl")
include("simulations.jl")

end # module
