using EqFlux
using Documenter
using Literate

cd(joinpath(@__DIR__, "src", "examples", "readme")) do
    Literate.markdown("readme.jl")
    Literate.notebook("readme.jl")
end

DocMeta.setdocmeta!(EqFlux, :DocTestSetup, :(using EqFlux); recursive=true)

makedocs(;
    modules=[EqFlux],
    authors="Ari Rappaport <ari.rappaport@inria.fr>",
    repo="https://github.com/triscale-innov/EqFlux.jl/blob/{commit}{path}#{line}",
    sitename="EqFlux.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aerappa.github.io/EqquilibratedFlux.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "examples/readme/readme.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/aerappa/EquilibratedFlux.jl",
    devbranch="main",
)
