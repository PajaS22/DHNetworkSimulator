using Documenter, DHNetworkSimulator

makedocs(sitename="DHNetworkSimulator", remotes = nothing)

deploydocs(
    repo = "github.com/PajaS22/DHNetworkSimulator.jl.git",
)
