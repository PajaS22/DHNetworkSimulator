# -----------------------------
# CONSTRUCTORS FOR Network
# -----------------------------

# create a network from a SimpleDiGraph
function Network(graph::SimpleDiGraph)
    # check that the graph is a tree
    if is_cyclic(graph) || !is_connected(graph)
        error("The input graph must be a tree (connected and acyclic).")
    end
    # check if there is at most one source node (node with outdegree > 0 and indegree == 0)
    source_nodes = [v for v in vertices(graph) if Graphs.outdegree(graph, v) > 0 && Graphs.indegree(graph, v) == 0]
    if length(source_nodes) > 1
        error("The input graph must have at most one source node (node with outdegree > 0 and indegree == 0). Found $(length(source_nodes)) source nodes.")
    end

    mg = MetaGraph(
        DiGraph();  # underlying graph structure
        label_type=String,
        vertex_data_type=NodeType,
        edge_data_type=EdgeType
    )
    default_labels = string.(1:nv(graph))
    for v in vertices(graph)
        mg[default_labels[v]] = EmptyNode()
    end
    for e in edges(graph)
        mg[default_labels[src(e)], default_labels[dst(e)]] = EmptyEdge()
    end

    return Network(mg, nothing, Set{String}(), Set{String}(), NeighborDicts())
end

"""Replace `EmptyEdge` placeholders with real `InsulatedPipe` objects.

`pipe_params` is a `Dict{Tuple{Int,Int}, PipeParams}` where each key is a
`(src_index, dst_index)` pair using the default integer labels assigned by
`Network(g::SimpleDiGraph)` (i.e. the vertex numbers from the original graph).
"""
function fill_physical_params!(network::Network, pipe_params::Dict{Tuple{Int,Int}, PipeParams})
    for (u, v) in keys(pipe_params)
        if has_edge(network, u, v)
            network[u, v] = InsulatedPipe(pipe_params[(u, v)])
        else
            error("Edge ($(u), $(v)) not found in the network graph.")
        end
    end
end

"""Set `(x, y)` coordinates for nodes in the network.

`node_positions` maps label → `(x, y)` tuple. Only nodes that are not `EmptyNode` are updated.
Positions are used by the visualization (`visualize_graph!`) but have no effect on simulation.
"""
function fill_node_positions!(network::Network, node_positions::Dict{String, Tuple{Float64, Float64}})
    for (label, pos) in node_positions
        if has_label(network, label)
            node_data = network[label]
            if !(node_data isa EmptyNode)
                node_data.common.position = pos
                network[label] = node_data  # update the node data in the network
            end
        else
            error("Node with label $(label) not found in the network.")
        end
    end
end

"""Rename the nodes in the network according to the provided `labels` vector.

The order of labels must match the order in which nodes were added to the graph.
When building from a `SimpleDiGraph`, the default labels are `"1"`, `"2"`, ..., `"n"`, so
`labels[1]` replaces `"1"`, `labels[2]` replaces `"2"`, and so on.
"""
function name_nodes!(network::Network, labels::Vector{String})
    # careful! The internal graph cannot rename the labels, so we have to remove a node and add it back again with the new label
    # tricky part is the removal, if we remove a node v, all the nodes v+k are reindexed so that the total |V| is reduced by 1
    # and the indexing doesnt contain gaps. The only consistent way to keep track of the nodes is by labels, they stay the same during removal
    # we use default names in form of "1", "2", ... "n" for the initial graph, so we can use the labels to find the correct node to rename
    
    if length(labels) != nv(network)
        error("Number of labels must match the number of nodes in the network.")
    end
    for (v, label) in zip(vertices(network), labels)
        rename_node!(network, string(v), label)
    end
end

"""Classify nodes as producer, loads, or junctions based on their connectivity.

Useful when building a network from a bare graph (e.g. loaded from a file) where node
types haven't been assigned yet. The classification rules are:

- `indegree == 0`, `outdegree > 0` → `ProducerNode` (the root / heat source)
- `indegree > 0`,  `outdegree == 0` → `LoadNode` (a leaf / consumer)
- `indegree > 0`,  `outdegree > 0` → `JunctionNode` (branching point)
- `indegree == 0`, `outdegree == 0` → `EmptyNode` (isolated — shouldn't appear in a valid network)

Also updates `network.producer_label` and `network.load_labels` accordingly.
"""
function identify_producer_and_loads!(network::Network)
    # identify producer node (source) and load nodes (sinks)
    for nl in all_labels(network)
        if outdegree(network, nl) > 0 && indegree(network, nl) == 0
            # this is a source node, we will treat it as producer
            if !(network[nl] isa ProducerNode)
                network[nl] = ProducerNode(nl)     
            end
        elseif outdegree(network, nl) == 0 && indegree(network, nl) > 0
            # this is a sink node, we will treat it as load
            if !(network[nl] isa LoadNode)
                network[nl] = LoadNode(nl)
            end
        elseif outdegree(network, nl) > 0 && indegree(network, nl) > 0
            # this is a junction node
            if !(network[nl] isa JunctionNode)
                network[nl] = JunctionNode(nl)
            end
        else
            # this is an isolated node, we will treat it as empty node
            warning("Node $(nl) is isolated (no incoming or outgoing edges). It will be treated as an empty node.")
            network[nl] = EmptyNode()
        end
    end
end

"""Convert nodes to `SumpNode` by label.

For each label in `sump_labels`, the node is converted to a `SumpNode` if it is currently
an `EmptyNode` or a `JunctionNode`. The node's `info` and `position` are preserved when
converting from a `JunctionNode`; a fresh `SumpNode` with `info=label` is created when
converting from an `EmptyNode`.

Nodes that are already a `SumpNode` are silently skipped. An error is raised if the label
does not exist in the network, or if the node is of a type that cannot be converted
(`ProducerNode`, `LoadNode`).

```julia
identify_sumps!(network, ["S1", "S2"])
identify_sumps!(network, "S1")          # single label convenience form
```
"""
function identify_sumps!(network::Network, sump_labels::Vector{String})
    for label in sump_labels
        has_label(network, label) || error("Node with label \"$label\" not found in the network.")
        node = network[label]
        if node isa SumpNode
            continue  # already a sump, nothing to do
        elseif node isa JunctionNode
            network[label] = SumpNode(NodeCommon(node.common.info, node.common.position, node.common.mass_flow))
        elseif node isa EmptyNode
            network[label] = SumpNode(label)
        else
            error("Cannot convert node \"$label\" ($(typeof(node))) to SumpNode. Only EmptyNode and JunctionNode can be converted.")
        end
    end
end
identify_sumps!(network::Network, label::String) = identify_sumps!(network, [label])

"""Set load power curve and relative mass-flow coefficient for each load node.

- `pwr_coefs`: `Dict{String, NTuple{3,Float64}}` mapping label → `(p₀, p₁, p₂)` for the quadratic
  demand curve `P(Tₐ) = p₀ + p₁·Tₐ + p₂·Tₐ²` (power in kW, temperature in °C).
- `m_r`: `Dict{String, Float64}` mapping label → relative mass-flow coefficient.

Every load node in the network must have an entry in both dictionaries.
The power curve is validated before storing — see [`validate_load_spec`](@ref).
"""
function fill_load_specs!(network::Network, pwr_coefs::Dict{String, NTuple{3, Float64}}, m_r::Dict{String, Float64})
    for label in network.load_labels
        haskey(pwr_coefs, label) || error("Power specifications for node $(label) not found in pwr_coefs.")
        haskey(m_r, label)       || error("Relative mass flow coefficient for node $(label) not found in m_r.")
        set_load_fn!(network, label, polynomial_load, collect(pwr_coefs[label]))
        set_load_m_rel!(network, label, m_r[label])
    end
end

"""Set the same load power curve and mass-flow coefficient for every load node in the network.

Handy for quick setups where all loads share the same parameters.
Defaults to zero power curve and `m_r=1.0` (equal flow split).
The power curve is validated before storing — see [`validate_load_spec`](@ref).
"""
function fill_load_specs!(network::Network; pwr_coefs = (0.0, 0.0, 0.0), m_r = 1.0)
    for label in network.load_labels
        set_load_fn!(network, label, polynomial_load, collect(pwr_coefs))
        set_load_m_rel!(network, label, m_r)
    end
end

"""Set the same load demand function and parameters for every load node in the network.

Use this overload when you need a function other than `polynomial_load`, e.g. `hockey_load`.

```julia
fill_load_specs!(network, hockey_load, [0.0, 0.05/15.0, 15.0]; m_r=1.0)
```
"""
function fill_load_specs!(network::Network, fn::Function, params::Vector{Float64}; m_r::Real = 1.0)
    for label in network.load_labels
        set_load_fn!(network, label, fn, params)
        set_load_m_rel!(network, label, Float64(m_r))
    end
end


function fill_load_specs!(network::Network, load::LoadSpec; m_r::Real = 1.0)
    for label in network.load_labels
        set_load_fn!(network, label, load.fn, load.params;
                     use_mass_flow=load.use_mass_flow, use_time=load.use_time)
        set_load_m_rel!(network, label, Float64(m_r))
    end
end