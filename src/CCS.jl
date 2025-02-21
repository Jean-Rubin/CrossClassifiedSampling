module CCS

using Distributions
using StatsBase

export PopParam, BaseModel, RatioPoisson, RatioIndep, RatioMixPoisson, RatioMixSimple
export PopParam3D, Model3D
export generate, simulate_completely, simulate_completely_var, simulate_completely_3d, simulate_completely_var_3d


include("models.jl")
include("sampling.jl")
include("bootstrap.jl")
include("intervals.jl")
include("simulations.jl")

end # module
