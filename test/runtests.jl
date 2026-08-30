using EquilibratedFlux
using Test

@testset "EquilibratedFlux.jl" begin
    include("dirichlettest.jl")
    include("neumanntest.jl")
    include("alloctest.jl")
    # Write your tests here.
end
