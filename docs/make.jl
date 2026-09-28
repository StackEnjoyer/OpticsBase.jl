using OpticsBase
using Documenter

DocMeta.setdocmeta!(OpticsBase, :DocTestSetup, :(using OpticsBase); recursive=true)

makedocs(;
    modules=[OpticsBase],
    authors="Hugo Uittenbosch <hugo.uittenbosch@dlr.de> and contributors",
    repo=Remotes.GitHub("StackEnjoyer", "OpticsBase.jl"),
    sitename="OpticsBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://StackEnjoyer.github.io/OpticsBase.jl",
        edit_link="main",
        assets=String[],
    ),
    pagesonly=true,
    pages=[
        "Home"             => "index.md",
        "Conventions"      => "conventions.md",
        "Ports"            => "ports.md",
        "Exchange formats" => "formats.md",
        "Solver interface" => "interface.md",
        "Index"            => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/StackEnjoyer/OpticsBase.jl.git",
    devbranch="main",
    push_preview=false,
)
