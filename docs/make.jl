using OrnsteinZernike
using Documenter

DocMeta.setdocmeta!(OrnsteinZernike, :DocTestSetup, :(using OrnsteinZernike); recursive=true)

makedocs(;
    modules=[OrnsteinZernike],
    authors="Edwin Bedolla",
    repo="https://github.com/edwinb-ai/OrnsteinZernike.jl/blob/{commit}{path}#{line}",
    sitename="OrnsteinZernike.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://edwinb-ai.github.io/OrnsteinZernike.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/edwinb-ai/OrnsteinZernike.jl",
    devbranch="main",
)
