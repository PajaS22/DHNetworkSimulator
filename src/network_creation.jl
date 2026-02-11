# -----------------------------
# CONSTRUCTORS FOR dhNetwork
# -----------------------------

# create a network from a SimpleDiGraph
function dhNetwork(graph::SimpleDiGraph)
    # check that the graph is a tree
    if is_cyclic(graph) || !is_connected(graph)
        error("The input graph must be a tree (connected and acyclic).")
    end
    # check if there is at most one source node (node with outdegree > 0 and indegree == 0)
    source_nodes = [v for v in vertices(graph) if outdegree(graph, v) > 0 && indegree(graph, v) == 0]
    if length(source_nodes) > 1
        error("The input graph must have at most one source node (node with outdegree > 0 and indegree == 0). Found $(length(source_nodes)) source nodes.")
    end

    mg = MetaGraph(
        DiGraph();  # underlying graph structure
        label_type=String,
        vertex_data_type=dhNodeType,
        edge_data_type=dhEdgeType
    )
    default_labels = string.(1:nv(graph))
    for v in vertices(graph)
        mg[default_labels[v]] = EmptyNode()
    end
    for e in edges(graph)
        mg[default_labels[src(e)], default_labels[dst(e)]] = EmptyEdge()
    end
    
    return dhNetwork{Int}(mg, nothing, Set{String}())
end

function fill_physical_params!(network::dhNetwork, pipe_params::Dict{Tuple{Int,Int}, PipeParams})
    for (u, v) in keys(pipe_params)
        if has_edge(network, u, v)
            network[u, v] = InsulatedPipe(pipe_params[(u, v)])
        else
            error("Edge ($(u), $(v)) not found in the network graph.")
        end
    end
end

function fill_node_positions!(network::dhNetwork, node_positions::Dict{String, Tuple{Float64, Float64}})
    for (label, pos) in node_positions
        if has_label(network, label)
            node_data = network[label]
            if !(node_data isa dhEmptyNode)
                node_data.common.position = pos
                network[label] = node_data  # update the node data in the network
            end
        else
            error("Node with label $(label) not found in the network.")
        end
    end
end

function name_nodes!(network::dhNetwork, labels::Vector{String})
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

function identify_producer_and_loads!(network::dhNetwork)
    # identify producer node (source) and load nodes (sinks)
    for nl in all_labels(network)
        if outdegree(network, nl) > 0 && indegree(network, nl) == 0
            # this is a source node, we will treat it as producer
            network[nl] = ProducerNode(nl)
        elseif outdegree(network, nl) == 0 && indegree(network, nl) > 0
            # this is a sink node, we will treat it as load
            network[nl] = LoadNode(nl)
        elseif outdegree(network, nl) > 0 && indegree(network, nl) > 0
            # this is a junction node
            network[nl] = JunctionNode(nl)
        else
            # this is an isolated node, we will treat it as empty node
            warning("Node $(nl) is isolated (no incoming or outgoing edges). It will be treated as an empty node.")
            network[nl] = EmptyNode()
        end
    end
end

function fill_load_specs!(network::dhNetwork, pwr_coefs::Dict{String, NTuple{3, Float64}}, m_r::Dict{String, Float64})
    for label in network.load_labels
        if haskey(pwr_coefs, label)
            set_load_pwr_coefs!(network, label, pwr_coefs[label])
        else
            error("Power specifications for node $(label) not found in load_coefs.")
        end
        if haskey(m_r, label)
            set_load_m_rel!(network, label, m_r[label])
        else
            error("Relative mass flow coefficient for node $(label) not found in m_r.")
        end
    end
end