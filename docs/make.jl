using EqFlux
using Documenter
using Literate

for name in ("readme", "Lshaped")
    cd(joinpath(@__DIR__, "src", "examples", name)) do
        Literate.markdown("$name.jl")
        Literate.notebook("$name.jl")
    end
end

DocMeta.setdocmeta!(EqFlux, :DocTestSetup, :(using EqFlux); recursive=true)

on_CI = get(ENV, "CI", "false") == "true"
makedocs(;
         modules=[EqFlux],
         authors="Ari Rappaport <ari.rappaport@inria.fr>",
         repo="https://github.comaerappa/EquilibratedFlux.jl/blob/{commit}{path}#{line}",
         sitename="EqFlux.jl",
         format=Documenter.HTML(;
                                prettyurls=get(ENV, "CI", "false") == "true",
                                canonical="https://aerappa.github.io/EqquilibratedFlux.jl",
                                edit_link="main",
                                assets=String[],
                                ),
         pages=[
             "Home" => "index.md",
             "Tutorials" => [
                 "examples/readme/readme.md",
                 "examples/Lshaped/Lshaped.md",
             ],
         ],
         warnonly = on_CI ? false : Documenter.except(:linkcheck_remotes),
         )

deploydocs(;
           repo="github.com/aerappa/EquilibratedFlux.jl",
           devbranch="main",
           )
