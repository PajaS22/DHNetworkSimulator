# API reference

This page collects the public API of DHNetworkSimulator in one place.

## Core types

### Network and topology

```@docs
DHNetworkSimulator.Network
DHNetworkSimulator.NeighborDicts
```

### Nodes and edges

```@docs
DHNetworkSimulator.NodeType
DHNetworkSimulator.EdgeType
DHNetworkSimulator.NodeCommon

DHNetworkSimulator.ProducerNode
DHNetworkSimulator.JunctionNode
DHNetworkSimulator.LoadNode
DHNetworkSimulator.EmptyNode

DHNetworkSimulator.PipeParams
DHNetworkSimulator.InsulatedPipe
DHNetworkSimulator.Plug
DHNetworkSimulator.EmptyEdge
```

### Simulation I/O

```@docs
DHNetworkSimulator.ProducerOutput
DHNetworkSimulator.SimulationResults
```

## Network construction and manipulation

```@docs
DHNetworkSimulator.fill_physical_params!
DHNetworkSimulator.fill_node_positions!
DHNetworkSimulator.fill_load_specs!
DHNetworkSimulator.name_nodes!
DHNetworkSimulator.identify_producer_and_loads!
```

```@docs
DHNetworkSimulator.has_label
DHNetworkSimulator.index_for
DHNetworkSimulator.rem_node!
DHNetworkSimulator.rename_node!
DHNetworkSimulator.vertices_data
DHNetworkSimulator.edges_data
DHNetworkSimulator.all_labels
```

## Graph interface

`Network` implements parts of the Graphs.jl interface (e.g., `nv`, `ne`, `vertices`, `edges`, `outdegree`, `indegree`, `neighbors`).

These label-based helpers are provided by DHNetworkSimulator:

```@docs
DHNetworkSimulator.outneighbors
DHNetworkSimulator.inneighbors
```

## Visualization and printing

```@docs
DHNetworkSimulator.visualize_graph!
DHNetworkSimulator.edge_info
DHNetworkSimulator.edge_info_hover
DHNetworkSimulator.edge_infos
DHNetworkSimulator.print_nodes
DHNetworkSimulator.print_edges
```

## Pipe geometry and hydraulics helpers

```@docs
DHNetworkSimulator.pipe_length
Base.length(::DHNetworkSimulator.InsulatedPipe)
DHNetworkSimulator.inner_diameter
DHNetworkSimulator.heat_resistance_forward
DHNetworkSimulator.heat_resistance_backward
DHNetworkSimulator.water_velocity
DHNetworkSimulator.water_velocities
```

## Plug-method helpers

```@docs
DHNetworkSimulator.collect_exiting_water_plugs!
DHNetworkSimulator.combine_plugs
DHNetworkSimulator.merge_same_temperature_plugs!
DHNetworkSimulator.merge_water_plug_vectors!
```

## Simulation

```@docs
DHNetworkSimulator.check_network!
DHNetworkSimulator.set_relative_mass_flows!
DHNetworkSimulator.set_absolute_mass_flows!
DHNetworkSimulator.steady_state_hydronynamics!
DHNetworkSimulator.fill_pipes_with_initial_temperature!
DHNetworkSimulator.time_step_thermal_dynamics!
DHNetworkSimulator.set_load_params!
DHNetworkSimulator.power_consumption
DHNetworkSimulator.consume_power!
DHNetworkSimulator.run_simulation
```

## Plotting

```@docs
DHNetworkSimulator.plot_simulation_results
```

## External helpers
This package relies on MetaGraphsNext.jl internally. If you need to convert a Graphs.jl vertex index back to the stored string label, use `MetaGraphsNext.label_for(nw.mg, idx)`.