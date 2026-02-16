# Distributed Heating Network Simulator

This is the documentation for **DHNetworkSimulator** ([source code on github](https://github.com/PajaS22/DHNetworkSimulator)).

DHNetworkSimulator is a Julia package for building, visualizing, and simulating **district heating networks**. A network is modeled as a **directed, acyclic tree**: a single producer sits at the root, water flows through junctions, and consumers (loads) are typically leaves.

The simulator uses a **quasi-dynamic** formulation to keep the problem tractable. At each time step, hydraulics (mass-flow distribution and pressures) are assumed to reach steady state much faster than temperature changes. Thermal dynamics are then dominated by advection (hot water moving through pipes) and heat losses to the ambient.

Each time step is therefore solved in two stages:
1. compute steady-state hydraulics,
2. advance temperatures with a **plug-flow** model that transports discrete “plugs” (parcels) of water through pipes. Plugs do not mix inside a pipe; mixing happens at junctions when flows merge and converge.

Simulations are typically **policy-driven**: you provide a `policy(t, Tₐ, T_back)` callback that returns the producer setpoints, and the simulator logs time series of temperatures, flows, and powers.

More details on the thermal model are in [Plug method](@ref).


## Getting started
To install the package, open the Julia REPL and type
```julia
julia> using Pkg; Pkg.add(url="https://github.com/PajaS22/DHNetworkSimulator")
```

### Examples

The easiest way to get a feel for the API is to run one of the example scripts. You can find them after cloning the [repo](https://github.com/PajaS22/DHNetworkSimulator) in subfolder scripts.
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

More examples can be found in [Simulation examples](@ref simulation_examples).