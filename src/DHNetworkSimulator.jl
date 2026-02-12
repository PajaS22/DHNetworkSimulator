module DHNetworkSimulator

using Graphs, MetaGraphsNext, DataStructures
using ForwardMethods                            # for overloading functions from Graphs.jl and MetaGraphs.jl on Network
using GraphMakie, Colors, ColorSchemes
using Plots
using Dates

if get(ENV, "GITHUB_ACTIONS", "false") != "true"
    using GLMakie # only load GLMakie when not running in GitHub Actions, to avoid issues with headless rendering
end

const WATER_DENSITY = 1000.0            # kg/m3 at approx. 25 degC
const WATER_SPECIFIC_HEAT = 4186.0      # J/(kg*K) at approx. 25 degC
const MINIMAL_RETURN_TEMPERATURE = 25.0 # °C, minimal return temperature for the simulation to avoid un-physical results

include("types.jl")

export NodeType, EdgeType, NodeCommon
export JunctionNode, LoadNode, ProducerNode, EmptyNode
export InsulatedPipe, EmptyEdge, Plug, PipeParams
export Network, NeighborDicts
export EmptyNode, EmptyEdge
export InsulatedPipe, JunctionNode, LoadNode, ProducerNode

include("network_creation.jl")
export fill_physical_params!
export name_nodes!
export identify_producer_and_loads!
export fill_node_positions!, fill_load_specs!

include("network_methods.jl")
export vertices_data, edges_data
export print_edges, print_nodes
export ne, nv, vertices, edges, all_labels
export index_for, has_label
export rem_node!, rename_node!
export length, inner_diameter, heat_resistance_forward, heat_resistance_backward
export water_velocities, water_velocity
export outneighbors, inneighbors, neighbors, degree, outdegree, indegree

include("printing.jl")

include("visualize_network.jl")
export visualize_graph!
export edge_info, edge_info_hover, edge_infos

include("simulation.jl")
export run_simulation
export set_relative_mass_flows!, set_absolute_mass_flows!
export steady_state_hydronynamics!
export fill_pipes_with_initial_temperature!
export time_step_thermal_dynamics!
export set_load_params!
export power_consumption, consume_power!

include("plot_simulation_results.jl")
export plot_simulation_results
export SimulationResults
export ProducerOutput

end # module DHNetworkSimulator
