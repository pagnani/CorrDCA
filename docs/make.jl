push!(LOAD_PATH,"../src/")
using CorrDCA
using Documenter
#DocMeta.setdocmeta!(CorrDCA, :DocTestSetup, :(using CorrDCA); recursive=true)
makedocs(;
    modules=[CorrDCA],
    authors="Andrea Pagnani",
    repo="https://github.com/pagnani/CorrDCA.jl/blob/{commit}{path}#L{line}",
    sitename="CorrDCA",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pagnani.github.io/CorrDCA",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pagnani/CorrDCA",
)
