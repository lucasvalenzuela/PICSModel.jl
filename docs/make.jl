using PICSModel
using Documenter

DocMeta.setdocmeta!(PICSModel, :DocTestSetup, :(using PICSModel); recursive=true)

makedocs(;
    modules=[PICSModel],
    authors="Lucas Valenzuela",
    sitename="PICSModel.jl",
    format=Documenter.HTML(;
        canonical="https://lucasvalenzuela.github.io/PICSModel.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/lucasvalenzuela/PICSModel.jl",
    devbranch="main",
)
