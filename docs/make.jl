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
        "Network Diffusion" => "network_diffusion.md",
        "SCC Tree" => "scctree.md",
        "Source → Sink Flows" => "source_sink.md",
        "Statistics" => "stats.md",
        "Export" => "export.md",
        "Utilities" => "utils.md",
    ],
)

deploydocs(
    repo = "github.com/alyst/HierarchicalHotNet.jl.git",
)
