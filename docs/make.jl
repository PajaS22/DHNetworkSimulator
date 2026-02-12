using Documenter, DHNetworkSimulator

makedocs(sitename="DHNetworkSimulator", remotes = nothing,
        workdir = joinpath(@__DIR__, ".."))

deploydocs(
    repo = "github.com/PajaS22/DHNetworkSimulator.jl.git",
)
