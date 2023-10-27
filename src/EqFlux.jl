module EqFlux
using Gridap
using Gridap.Geometry
#using BenchmarkTools

include("SpaceManager.jl")
include("CellwiseAssembler.jl")
include("Patch.jl")
include("Util.jl")
include("DOFManager.jl")
include("LinAlgAssembler.jl")
include("AveragedFlux.jl")
include("FluxBuilder.jl")

export build_equilibrated_flux
export build_averaged_flux


end # module
