using CorrDCA
using Documenter

makedocs(;
    modules=[CorrDCA],
    authors="Andrea Pagnani",
    repo="https://github.com/pagnani/CorrDCA.jl/blob/{commit}{path}#L{line}",
    sitename="CorrDCA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pagnani.github.io/CorrDCA.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pagnani/CorrDCA.jl",
)
