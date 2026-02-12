# DHNetworkSimulator

DHNetworkSimulator is a Julia package for building, visualizing, and simulating **district heating networks**. Network is represented as directed tree graph, root node is producer and leaf are consumers (loads).

Thermodynamics is solved using quazi-dynamic assumption, which makes the problem tractable. The assumption is that hydraulics (change in flow, pressures) has frequency much higher than thermics (change in temperature). In simple words, pressure moves in water with speed of sound and temperature moves mainly with advection (the hot water itself moves).

We solve the dynamics in two steps:
1. find steady state of hydraulics
2. solve thermics by moving "plugs" of water. Consecutive plugs dont mix.

More on plug method in [Plug method](Plug_method.md).

## Documentation

See documentation at (https://pajas22.github.io/DHNetworkSimulator)

## Features

- Network model built on `Graphs.jl` + `MetaGraphsNext.jl`
- Node/edge types for producers, loads, junctions, and pipes
- Steady-state hydrodynamics solver (`steady_state_hydronynamics!`)
- Time-stepping thermal simulation (`run_simulation`)
- feedback control on the producer side (using policy function)
- Visualization (`visualize_graph!`) and result plotting (`plot_simulation_results`)

## Requirements

- Julia 1.11+ (tested with Julia 1.12)

## Installation

Clone the repository and instantiate the project environment:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Then load the package:

```julia
using DHNetworkSimulator
```

## Quick start

The easiest way to get a feel for the API is to run one of the example scripts.
- `scripts/basic_network.jl`
- `scripts/bigger_network.jl`



### Minimal example (policy-driven simulation)

`run_simulation` takes a network, a time vector, and a `policy(t, Tₐ, T_back)` callback that returns a `ProducerOutput`. Simulation than returns `SimulationResults` struct holding  time series of many physical variables of the network.

```julia
using DHNetworkSimulator

network = Network()
network["producer"] = ProducerNode((0.0, 0.0))

# ... add junctions/loads/pipes here ...

t = collect(0.0:60.0:3*3600.0) # seconds

function policy(t, Tₐ, T_back) # constant flow and temperature
	return ProducerOutput(mass_flow=15.0, temperature=90.0)
end

# set initial temperatures in pipes to 75°C and 60°C
results = run_simulation(network, t, policy; T0_f=75.0, T0_b=60.0)
# plot temperature of water entering loads
plot_simulation_results(results, :T_load_in)
```

## Data

- `datasets/` contains sample network definition with physical parameters of the pipes and load characteristics of the consumers. There are also real measurements of outdoor temperature for duration of one year.
- this network can be simulated as shown in script `scripts/bigger_network.jl` 


## Project structure

- `src/` — package source code
- `scripts/` — runnable examples and experiments
- `datasets/` — sample network + outdoor temperature
- `test/` — unit tests

## License

See `LICENSE`.
