using DistributionMetrics
using Documenter

DocMeta.setdocmeta!(DistributionMetrics, :DocTestSetup, :(using DistributionMetrics); recursive=true)

makedocs(;
    modules=[DistributionMetrics],
    authors="Šimon Soldát",
    sitename="DistributionMetrics.jl",
    format=Documenter.HTML(;
        canonical="https://soldasim.github.io/DistributionMetrics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/soldasim/DistributionMetrics.jl",
    devbranch="main",
)
