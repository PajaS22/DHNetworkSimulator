using Test
using DHNetworkSimulator
using MetaGraphsNext, Graphs

const DHN = DHNetworkSimulator

include("network_tests.jl")
include("constructor_tests.jl")
include("network_creation_tests.jl")
include("physics_tests.jl")
include("simulation_tests.jl")