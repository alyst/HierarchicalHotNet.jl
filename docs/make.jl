using Revise
using HierarchicalHotNet
using Documenter

cd("docs")

makedocs(
    root = pwd(),#joinpath(dirname(pathof(HierarchicalHotNet)), "..", "docs"),
    format = Documenter.HTML(prettyurls=false),
    sitename = "HierarchicalHotNet.jl",
    authors = "Alexey Stukalov",
    modules = [HierarchicalHotNet],
    clean = true,
    pages = [
        "Introduction" => "index.md",
        "SCC Tree" => "scctree.md",
    ],
)

deploydocs(
    repo = "github.com/alyst/HierarchicalHotNet.jl.git",
)
