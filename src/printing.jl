# ------------------------------------------------- #
# NICE PRINTING METHODS FOR DIFFERENT TYPES
# ------------------------------------------------- #

# ------------------------------------------------ #
# nice printing of Network
# ------------------------------------------------ #

function Base.show(io::IO, nw::Network)
    println(io, "DH Network:")
    println(io, " Number of nodes: ", nv(nw))
    println(io, " Number of edges: ", ne(nw))
    if !ismissing(nw.producer_label)
        println(io, " Producer node: ", nw.producer_label)
    else
        println(io, " Producer node: not set")
    end
    if !ismissing(nw.load_labels)
        println(io, " Load nodes ($(length(nw.load_labels))): {", join(nw.load_labels, ", "), "}")
    else
        println(io, " Load nodes: not set")
    end
    # print nodes
    println(io, "---")
    print_nodes(nw)
    # print edges
    println(io, "---")
    print_edges(nw)
end
# print edges(network)
function print_edges(nw::Network)
    println("Edges in the network:")
    for e in edges(nw)
        src_label = label_for(nw, src(e))
        dst_label = label_for(nw, dst(e))
        print("[$(src(e))] $src_label -->  [$(dst(e))] $dst_label : ", nw[src_label, dst_label])
    end
end
# print nodes(network)
function print_nodes(nw::Network)
    println("Nodes in the network:")
    for v in labels(nw.mg)
        print("[$(code_for(nw.mg, v))] $v : ", nw.mg[v])
    end
end

# ------------------------------------------------ #
# nice printing of Node
# ------------------------------------------------ #
function string_node_common(common::NodeCommon)::String
    summary = "Info: $(common.info)"
    summary *= (!ismissing(common.mass_flow) ? ", Mass flow: $(round(common.mass_flow; digits=2)) kg/s" : "")
    return summary
end
function Base.show(io::IO, node::NodeType)
    print(io, "$(typeof(node))")
end
# Empty node (minimal)
function Base.show(io::IO, ::EmptyNode)
    println(io, "Empty Node")
end

# Junction node (with details)
function Base.show(io::IO, node::JunctionNode)
    summary = "Junction Node, "
    summary *= string_node_common(node.common)
    println(io, summary)
end

# Load node (with load info)
function Base.show(io::IO, node::LoadNode)
    summary = "Load Node, "
    summary *= string_node_common(node.common)
    summary *= !ismissing(node.load) ? ", Load (at 0°C): $(round(node.load[1]; digits=1)) kW" : ""
    summary *= !ismissing(node.m_rel) ? ", M_flow_rel: $(round(node.m_rel; digits=2))" : ""
    println(io, summary)
end

# Producer node (with details)
function Base.show(io::IO, node::ProducerNode)
    summary = "Producer Node, "
    summary *= string_node_common(node.common)
    println(io, summary)
end

# ------------------------------------------------ #
# nice printing of Edge
# ------------------------------------------------ #

function Base.show(io::IO, edge::EdgeType)
    print(io, "$(typeof(edge))")
end
# Empty edge (minimal)
function Base.show(io::IO, ::EmptyEdge)
    println(io, "Empty Edge")
end
# Pipe edge (with details)
function Base.show(io::IO, edge::InsulatedPipe)
    summary = "Pipe Edge, L=$(pipe_length(edge)), D_in=$(inner_diameter(edge)), R_f=$(heat_resistance_forward(edge)), R_b=$(heat_resistance_backward(edge))"
    summary *= !ismissing(edge.mass_flow) ? ", Mass flow: $(round(edge.mass_flow; digits=2)) kg/s" : ""
    summary *= !ismissing(edge.m_rel) ? ", M_flow_rel: $(round(edge.m_rel; digits=2))" : ""
    if (!isempty(edge.plugs_f))
        summary *= ", Plugs (forward): ["
        for (i, plug) in enumerate(edge.plugs_f)
            summary *= "(T=$(round(plug.T;digits=1)), m=$(round(plug.m; digits=2)))"
            if i < length(edge.plugs_f)
                summary *= ", "
            end
        end
        summary *= "]"
    end
    println(io, summary)
end
