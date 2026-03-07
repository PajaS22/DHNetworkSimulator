# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

**DHNetworkSimulator** is a Julia package for building, visualizing, and simulating district heating (DH) networks. It uses a quasi-dynamic approach: hydraulics (mass flow) solved at steady state, thermal dynamics via a plug-flow method where water is represented as discrete parcels.

## Commands

```bash
# Install dependencies
julia --project -e "using Pkg; Pkg.instantiate()"

# Run all tests
julia --project -e "using Pkg; Pkg.test()"

# Run a single test file
julia --project test/network_tests.jl
julia --project test/constructor_tests.jl
julia --project test/network_creation_tests.jl
julia --project test/physics_tests.jl
julia --project test/simulation_tests.jl

# Run example scripts
julia --project scripts/basic_network.jl
```

## Architecture

### Source Files (`src/`)

| File | Purpose |
|------|---------|
| `DHNetworkSimulator.jl` | Module root: imports, exports, constants (`WATER_DENSITY`, `WATER_SPECIFIC_HEAT`, `MINIMAL_RETURN_TEMPERATURE`) |
| `types.jl` | All data structures: `Network`, node types, edge types, `Plug`, `PipeParams` |
| `network_creation.jl` | Incremental network construction helpers (`fill_physical_params!`, `identify_producer_and_loads!`, etc.) |
| `network_methods.jl` | Label-based indexing, neighbor queries, topology modification, `check_network!` |
| `simulation.jl` | Hydraulics (`steady_state_hydronynamics!`) and thermal simulation (`run_simulation`) |
| `visualize_network.jl` | GraphMakie interactive visualization (`visualize_graph!`) |
| `plot_simulation_results.jl` | Plots.jl time-series plots for `SimulationResults` |
| `printing.jl` | `Base.show` overloads for network types |

### Key Design Decisions

**Label-based indexing**: Nodes and edges are always accessed by string labels, not integer indices. Integer indices are an internal MetaGraphsNext detail.

**Lazy neighbor cache**: `NeighborDicts` in `Network` is rebuilt on demand. Any topology modification must mark it dirty; `check_and_update_neighbor_dicts!` rebuilds it.

**Two plug queues per pipe**: Both `InsulatedPipe` and `ZeroPipe` have `plugs_f` (supply/forward) and `plugs_b` (return/backward) fields. These are `Vector{Plug}` — front of the vector is the pipe outlet. For `ZeroPipe` the queues act as within-timestep pass-through buffers and are always empty at the start/end of each step.

**Tree-only topology**: `check_network!` (called at simulation start) enforces that the graph is a DAG, has exactly one `ProducerNode`, and all nodes connect back to the producer.

**Two-phase hydraulics**: First compute `m_rel` (relative split coefficients) via post-order DFS, then compute absolute `mass_flow` via BFS from root.

**Load power model**: `LoadNode.load` holds a `LoadSpec(fn, params)` where `fn(params, T_a) -> kW`. The built-in `polynomial_load` gives the quadratic default. Custom functions set via `set_load_fn!`; parameters updated via `set_load_params!`. Both validate non-negativity and finiteness over a configurable T_a range before storing.

### Node Types

- `ProducerNode` — root of the tree; at most one per network
- `LoadNode` — leaf nodes (outdegree=0) with power demand curves
- `JunctionNode` — branching/merging points (indegree≥1, outdegree≥1)
- `EmptyNode` — placeholder used during construction

### Edge Types

- `InsulatedPipe` — physical pipe with `PipeParams` (length, inner diameter, thermal resistances) and two plug queues (`plugs_f`, `plugs_b`). All thermal/hydraulic accessors (`pipe_length`, `inner_diameter`, `heat_resistance_forward/backward`, `water_velocity`) operate on this type.
- `ZeroPipe` — a distinct `mutable struct <: EdgeType` modelling a zero-delay, lossless connection. Carries no `physical_params`; all dimension accessors return `0.0` and `water_velocity` returns `missing`. Use `e isa ZeroPipe` to identify it. Useful for the zero-delay branch of an aggregated junction pair.
- `EmptyEdge` — immutable placeholder used during topology construction; replace before simulating.

### Simulation Flow (per time step)

1. `steady_state_hydronynamics!` — recompute mass flows
2. `time_step_thermal_dynamics_forward!` — BFS from producer, advect supply plugs
3. Load power consumption cools plugs at leaves
4. `time_step_thermal_dynamics_backward!` — DFS from leaves, merge return flows at junctions
5. Heat losses applied exponentially along pipes

### Policy Function Signature

```julia
# Standard (forward + backward)
function policy(t, T_ambient, T_back) -> ProducerOutput(mass_flow=..., temperature=...)

# forward_only=true
function policy(t, T_ambient) -> ProducerOutput(mass_flow=..., temperature=...)
```

### GLMakie Environment Guard

`DHNetworkSimulator.jl` skips loading GLMakie when `GITHUB_ACTIONS=true` is set in the environment (for headless GitHub Actions). Visualization functions require GLMakie.
