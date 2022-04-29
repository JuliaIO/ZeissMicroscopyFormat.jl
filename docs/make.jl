using ZeissMicroscopyFormat
using Documenter

DocMeta.setdocmeta!(ZeissMicroscopyFormat, :DocTestSetup, :(using ZeissMicroscopyFormat); recursive=true)

makedocs(;
    modules=[ZeissMicroscopyFormat],
    authors="Tim Holy <tim.holy@gmail.com> and contributors",
    repo="https://github.com/JuliaIO/ZeissMicroscopyFormat.jl/blob/{commit}{path}#{line}",
    sitename="ZeissMicroscopyFormat.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaIO.github.io/ZeissMicroscopyFormat.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaIO/ZeissMicroscopyFormat.jl",
    devbranch="main",
)
