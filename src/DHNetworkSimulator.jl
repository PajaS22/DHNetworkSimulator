module DHNetworkSimulator

using Graphs, MetaGraphsNext, DataStructures
using ForwardMethods                            # for overloading functions from Graphs.jl and MetaGraphs.jl on Network
using GraphMakie, Colors, ColorSchemes
using Plots
using Dates
using PrecompileTools

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
export EmptyNode, EmptyEdge, ZeroPipe
export InsulatedPipe, JunctionNode, LoadNode, ProducerNode
export LoadSpec, polynomial_load, hockey_load, DEFAULT_LOAD_PARAMS

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
export rem_node!, rename_node!, remove_edge!
export pipe_length, length, inner_diameter, heat_resistance_forward, heat_resistance_backward, mass_flow, m_rel, volume
export water_velocities, water_velocity
export outneighbors, inneighbors, neighbors, degree, outdegree, indegree
export validate_load_spec, set_load_fn!


include("printing.jl")

include("visualize_network.jl")
export visualize_graph!
export NodeHighlight, highlight_nodes!, reset_highlights!
export edge_info, edge_info_hover, edge_infos
export compute_zero_pipe_load_positions
export DEFAULT_ZERO_PIPE_K_ATTRACTION, DEFAULT_ZERO_PIPE_K_REPULSION

include("simulation.jl")
export run_simulation
export set_relative_mass_flows!, set_absolute_mass_flows!
export steady_state_hydrodynamics!
export fill_pipes_with_initial_temperature!
export time_step_thermal_dynamics!
export set_load_params!
export power_consumption, consume_power!

include("plot_simulation_results.jl")
export plot_simulation_results
export SimulationResults
export ProducerOutput

include("network_analysis.jl")
export producer_loads_volumes, producer_loads_delays

# Precompile simulation and auto-positioning so the first user call pays no JIT cost.
@compile_workload begin
    # ── simulation ────────────────────────────────────────────────────────────
    let
        # 4-node network: producer → junction → 2 loads (mirrors the test fixture)
        nw = Network()
        nw["_wm_p"]  = ProducerNode((0.0, 0.0))
        nw["_wm_j"]  = JunctionNode((50.0, 0.0))
        nw["_wm_l1"] = LoadNode("_wm_l1", (100.0,  50.0), 1.0)
        nw["_wm_l2"] = LoadNode("_wm_l2", (100.0, -50.0), 2.0)
        nw["_wm_p",  "_wm_j"]  = InsulatedPipe(50.0)
        nw["_wm_j",  "_wm_l1"] = InsulatedPipe(50.0)
        nw["_wm_j",  "_wm_l2"] = InsulatedPipe(50.0)
        run_simulation(nw, [0.0, 60.0, 120.0],
            (t, Ta, Tb) -> ProducerOutput(mass_flow=5.0, temperature=80.0))
    end
    # ── ZeroPipe auto-positioning ─────────────────────────────────────────────
    let
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2); add_edge!(g, 2, 3)
        nw = Network(g)
        name_nodes!(nw, ["_wm_pp", "_wm_pj", "_wm_pl"])
        identify_producer_and_loads!(nw)
        fill_node_positions!(nw, Dict("_wm_pp" => (0.0, 0.0), "_wm_pj" => (100.0, 0.0)))
        nw["_wm_pp", "_wm_pj"] = InsulatedPipe(100.0)
        nw["_wm_pj", "_wm_pl"] = ZeroPipe()
        compute_zero_pipe_load_positions(nw.mg)
    end
end

# Precompile visualization — only when GLMakie is loaded (not in headless CI).
# Wrapped in try-catch so precompilation survives on headless machines that are
# not GitHub Actions (e.g. GPU-less servers, WSL without a display).
if get(ENV, "GITHUB_ACTIONS", "false") != "true"
    @compile_workload begin
        let
            nw = Network()
            nw["_wv_p"]  = ProducerNode((0.0, 0.0))
            nw["_wv_j"]  = JunctionNode((50.0, 0.0))
            nw["_wv_l1"] = LoadNode("_wv_l1", (100.0,  50.0), 1.0)
            nw["_wv_l2"] = LoadNode("_wv_l2", (100.0, -50.0), 2.0)
            nw["_wv_p",  "_wv_j"]  = InsulatedPipe(50.0)
            nw["_wv_j",  "_wv_l1"] = InsulatedPipe(50.0)
            nw["_wv_j",  "_wv_l2"] = InsulatedPipe(50.0)
            try
                visualize_graph!(nw)
            catch
            end
        end
    end
end

end # module DHNetworkSimulator
