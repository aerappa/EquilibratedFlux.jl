using EquilibratedFlux
using Documenter
using Literate

for name in ("readme", "Lshaped")
    cd(joinpath(@__DIR__, "src", "examples", name)) do
        Literate.markdown("$name.jl")
        Literate.notebook("$name.jl")
        if name == "readme"
            @info "running script `$name.jl`"
            include(joinpath(pwd(), "$name.jl"))
        end
    end
end

DocMeta.setdocmeta!(EquilibratedFlux, :DocTestSetup, :(using EquilibratedFlux); recursive=true)

on_CI = get(ENV, "CI", "false") == "true"
makedocs(;
         modules=[EquilibratedFlux],
         authors="Ari Rappaport <ari.rappaport@inria.fr>",
         repo="https://github.com/aerappa/EquilibratedFlux.jl/blob/{commit}{path}#{line}",
         sitename="EquilibratedFlux.jl",
         format=Documenter.HTML(;
                                prettyurls=get(ENV, "CI", "false") == "true",
                                canonical="https://aerappa.github.io/EquilibratedFlux.jl",
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
