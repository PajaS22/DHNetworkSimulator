using Documenter, DHNetworkSimulator

makedocs(sitename="DHNetworkSimulator", remotes = nothing,
        workdir = joinpath(@__DIR__, ".."),
        pages = [
        "Home" => "index.md",
        "Plug Method" => "Plug_method.md",
        "Network" => ["Network" => "Network.md", "Nodes and Edges" => "types.md"],
        "Simulation" => "simulation.md",
        ])

deploydocs(
    repo = "github.com/PajaS22/DHNetworkSimulator.jl.git",
)
