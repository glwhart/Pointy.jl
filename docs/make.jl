using Pointy
using Documenter

DocMeta.setdocmeta!(Pointy, :DocTestSetup, :(using Pointy); recursive=true)

makedocs(;
    modules=[Pointy],
    authors="Gus Hart",
    repo="https://github.com/glwhart/Pointy.jl/blob/{commit}{path}#{line}",
    sitename="Pointy.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://glwhart.github.io/Pointy.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/glwhart/Pointy.jl",
)
