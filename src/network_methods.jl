# ------------------------------------------------ #
# METHODS FOR NETWORK MANIPULATION AND ANALYSIS
# ------------------------------------------------ #
# we will overload some functions from Graphs.jl and MetaGraphs.jl for use on Network's MetaGraph

import Graphs: nv, ne, vertices, edges

# ------------------------------------------------ #
# accessing network node and edge data
# ------------------------------------------------ #

# indexing nodes
Base.getindex(nw::Network, label::String)::NodeType = nw.mg[label]
Base.getindex(nw::Network, labels::Vector{String}) = [nw.mg[lbl] for lbl in labels]
Base.getindex(nw::Network, i::Int)::NodeType = nw.mg[MetaGraphsNext.label_for(nw.mg, i)]
Base.getindex(nw::Network, i::Vector{Int}) = [nw.mg[MetaGraphsNext.label_for(nw.mg, j)] for j in i]

Base.setindex!(nw::Network, v::ProducerNode, label::String) = add_producer_node!(nw, v, label)
Base.setindex!(nw::Network, v::LoadNode, label::String) = add_load_node!(nw, v, label)
Base.setindex!(nw::Network, v::NT, label::String) where {NT<:NodeType} = add_node!(nw, v, label)
Base.setindex!(nw::Network, v::Vector{NT}, labels::Vector{String}) where {NT<:NodeType} = (for (j, lbl) in enumerate(labels); nw.mg[lbl] = v[j]; end)

"""Return `true` if `label` exists as a node label in the network."""
has_label(nw::Network, label::String) = MetaGraphsNext.haskey(nw.mg, label)

"""Return the internal vertex index for a given node `label`.

This is mainly useful when interacting with Graphs.jl functions that operate on
integer vertex indices.
"""
index_for(nw::Network, label::String) = MetaGraphsNext.code_for(nw.mg, label)
@forward_methods Network field=mg MetaGraphsNext.label_for(_,i::Int)

# indexing edges
has_edge(nw::Network, src::Int, dst::Int) = Graphs.has_edge(nw.mg, src, dst)
# has_edge(nw::Network, src::Int, dst::Int) = has_edge(nw, label_for(nw.mg, src), label_for(nw.mg, dst))

Base.setindex!(nw::Network, e::ET, src::String, dst::String) where {ET<:EdgeType} = (nw.mg[src, dst] = e)
function index_for(nw::Network, src::String, dst::String)
    if(!Graphs.has_edge(nw.mg, src, dst))
        error("Edge does not exist between $src and $dst")
    end
    return (MetaGraphsNext.code_for(nw.mg, src), MetaGraphsNext.code_for(nw.mg, dst))
end
Base.getindex(nw::Network, src::String, dst::String)::EdgeType = nw.mg[src, dst]
function Base.getindex(nw::Network, v_i::Int, u_i::Int)
    k1 = label_for(nw.mg, v_i)
    k2 = label_for(nw.mg, u_i)
    if MetaGraphsNext.haskey(nw.mg, k1, k2)
        return nw.mg[k1, k2]
    else
        error("Edge does not exist between vertices $v_i and $u_i")
    end
end
Base.setindex!(nw::Network, e::ET, v_i::Int, u_i::Int)  where {ET<:EdgeType} = (nw.mg[label_for(nw.mg,v_i), label_for(nw.mg,u_i)] = e)

function check_and_update_neighbor_dicts!(nw::Network)
    if !nw.neighbor_dicts.need_rebuild # flag is off, no need to rebuild
        return false
    end
    for v in vertices(nw.mg)
        nw.neighbor_dicts.outneighbors[label_for(nw.mg, v)] = [label_for(nw.mg, w) for w in Graphs.outneighbors(nw.mg, v)]
        nw.neighbor_dicts.inneighbors[label_for(nw.mg, v)] = [label_for(nw.mg, w) for w in Graphs.inneighbors(nw.mg, v)]
    end
    nw.neighbor_dicts.need_rebuild = false # reset the flag after rebuilding
    return true
end

"""Check network upon start of the simulation.

1. update neighbor dicts if needed
2. if there was a change, check that
    - there is exactly one producer node
    - there are no cycles in the network (DAG)
    - all nodes are reachable from the producer node
    - load nodes are leaves
"""
function check_network!(network::Network)
    if check_and_update_neighbor_dicts!(network) # there was a change
        # check that there is exactly one producer node
        if isnothing(network.producer_label)
            error("Producer node is not set in the network!")
        end
        # check that there are no cycles in the network (DAG)
        if is_cyclic(network.mg)
            error("The network graph must be acyclic!")
        end
        # check that all nodes are reachable from the producer node
        if !is_connected(network.mg)
            error("The network graph must be connected!")
        end
        # check that load nodes are leaves
        for load_label in network.load_labels
            if outdegree(network, load_label) != 0
                error("Load nodes must be leaves (outdegree must be 0)! Node with label $load_label has outdegree $(outdegree(network, load_label))")
            end
        end
    end
end

# removing nodes and edges
"""Remove a node from a `Network` by its string label.

This updates `producer_label` / `load_labels` as needed and marks neighbor caches
dirty.
"""
function rem_node!(nw::Network, label::String)
    v = nw[label]
    if v isa ProducerNode
        # removing producer node
        nw.producer_label = nothing
    elseif v isa LoadNode
        # removing load node
        filter!(x -> x != label, nw.load_labels)
    end
    Graphs.rem_vertex!(nw.mg, index_for(nw, label))
    nw.neighbor_dicts.need_rebuild = true
end
rem_node!(nw::Network, idx::Integer) = rem_node!(nw, label_for(nw, idx))

# getting all node data and edge data
"""Return a vector of all node data stored in a `Network` (in label order)."""
vertices_data(nw::Network) = [nw.mg[v] for v in labels(nw.mg)]

"""Return a vector of all node data stored in a MetaGraphsNext `MetaGraph` (in label order)."""
vertices_data(mg::MetaGraph) = [mg[v] for v in labels(mg)]

"""Return a vector of all edge data stored in a `Network`."""
edges_data(nw::Network) = [nw.mg[src,dst] for (src,dst) in MetaGraphsNext.edge_labels(nw.mg)]

"""Return all vertex labels in the network."""
all_labels(nw::Network) = [l for l in MetaGraphsNext.labels(nw.mg)]


# forwarding Network indexing to its MetaGraph
@forward_methods Network field=mg nv ne vertices edges
@forward_methods Network field=mg Graphs.degree(_,i::Int) Graphs.outdegree(_,i::Int) Graphs.inneighbors(_,i::Int) Graphs.outneighbors(_,i::Int)

# my implementation of neighbors and degree functions using the neighbor dicts for efficient access during simulation
"""Return outgoing neighbor labels for a node label.

This is the `Network`-specific overload that returns **string labels** (not
Graphs.jl vertex indices). Neighbor caches are rebuilt automatically when
needed.
"""
function outneighbors(nw::Network, label::String)
    check_and_update_neighbor_dicts!(nw)
    nw.neighbor_dicts.outneighbors[label]
end
"""Return incoming neighbor labels for a node label.

This is the `Network`-specific overload that returns **string labels** (not
Graphs.jl vertex indices). Neighbor caches are rebuilt automatically when
needed.
"""
function inneighbors(nw::Network, label::String)
    check_and_update_neighbor_dicts!(nw)
    nw.neighbor_dicts.inneighbors[label]
end
degree(nw::Network, label::String)  = outdegree(nw, label) + indegree(nw, label)
outdegree(nw::Network, label::String) = length(outneighbors(nw, label))
indegree(nw::Network, label::String) = length(inneighbors(nw, label))


# Node data operations: positions
position(::EmptyNode) = missing
position(v::NT)  where {NT<:NodeType} = !ismissing(v.common.position) ? (v.common.position[1], v.common.position[2]) : missing
positions(mg::MetaGraph) = [position(v) for v in vertices_data(mg)]
distance(v1::NodeType, v2::NodeType) = sqrt((v1.common.position[1] - v2.common.position[1])^2 + (v1.common.position[2] - v2.common.position[2])^2)

# ------------------------------------------------ #
# NODE TYPE SPECIFIC METHODS
# ------------------------------------------------ #

function add_load_node!(nw::Network, v::LoadNode, label::String)
    # check if node with the same label already exists
    add_node!(nw, v, label)
    push!(nw.load_labels, label)
end
function add_producer_node!(nw::Network, v::ProducerNode, label::String)
    if(!isnothing(nw.producer_label) && nw.producer_label != label)
        error("Producer node is already set in the network! (label: $(nw.producer_label)) Remove it first or rewrite it.")
    end

    add_node!(nw, v, label)
    nw.producer_label = label
end

# default method for adding a node in network
function add_node!(nw::Network, v::NT, label::String) where {NT<:NodeType}
    if has_label(nw, label)
        if nw[label] isa LoadNode
            # removing load node from load_labels
            filter!(x -> x != label, nw.load_labels)
        elseif nw[label] isa ProducerNode
            # removing producer node
            nw.producer_label = nothing
        end
    end
    nw.mg[label] = v
    nw.neighbor_dicts.need_rebuild = true
end

"""Rename a node label in-place.

This preserves the node data and all incident edges by removing the old label and
re-inserting under `new_label`.
"""
function rename_node!(nw::Network, old_label::String, new_label::String)
    if !has_label(nw, old_label)
        error("Node with label $old_label does not exist in the network.")
    end
    if has_label(nw, new_label)
        error("Node with label $new_label already exists in the network.")
    end
    node = nw[old_label]
    in_edges = [(label, nw[label, old_label]) for label in DHNetworkSimulator.inneighbors(nw, old_label)]
    out_edges = [(label, nw[old_label, label]) for label in DHNetworkSimulator.outneighbors(nw, old_label)]
    rem_node!(nw, old_label)
    nw[new_label] = node
    for (src_label, edge_data) in in_edges
        nw[src_label, new_label] = edge_data
    end
    for (dst_label, edge_data) in out_edges
        nw[new_label, dst_label] = edge_data
    end
    nw.neighbor_dicts.need_rebuild = true
    check_and_update_neighbor_dicts!(nw) # update the neighbor dicts after renaming
end

# set input mass flow and temperature at source node
function set_source_inputs!(nw::Network, mass_flow::Float64, temperature::Float64)
    if isnothing(nw.producer_label)
        error("Source node is not set in the network!")
    end
    src_node = nw[nw.producer_label]
    src_node.common.mass_flow = mass_flow
    src_node.common.temperature = temperature
end

# set coefficients for load node
function set_load_pwr_coefs!(nw::Network, label::String, coefs::Tuple{Float64, Float64, Float64})
    if !has_label(nw, label)
        error("Node with label $label does not exist in the network.")
    end
    node = nw[label]
    if !(node isa LoadNode)
        error("Node with label $label is not a load node.")
    end
    node.load = coefs
    nw[label] = node
end

# set relative_mass_flow for load node
function set_load_m_rel!(nw::Network, label::String, m_rel::Float64)
    if !has_label(nw, label)
        error("Node with label $label does not exist in the network.")
    end
    node = nw[label]
    if !(node isa LoadNode)
        error("Node with label $label is not a load node.")
    end
    node.m_rel = m_rel
    nw[label] = node
end

# ------------------------------------------------ #
# EDGE TYPE SPECIFIC METHODS
# ------------------------------------------------ #

"""Return physical pipe length in meters."""
pipe_length(e::InsulatedPipe) = e.physical_params.length

"""Alias for [`pipe_length`](@ref).

This makes `length(pipe)` meaningful for `InsulatedPipe` and matches the package export list.
"""
Base.length(e::InsulatedPipe) = pipe_length(e)

"""Return pipe inner diameter in meters."""
inner_diameter(e::InsulatedPipe) = e.physical_params.inner_diameter

"""Return thermal resistance (supply direction) in m·K/W."""
heat_resistance_forward(e::InsulatedPipe) = e.physical_params.heat_resistance_forward

"""Return thermal resistance (return direction) in m·K/W."""
heat_resistance_backward(e::InsulatedPipe) = e.physical_params.heat_resistance_backward

"""Compute water velocity in an `InsulatedPipe` in m/s.

Returns `missing` if the pipe has no mass flow set.
"""
function water_velocity(e::InsulatedPipe)
    if(e.mass_flow === missing)
        return missing
    end
    if(inner_diameter(e) === missing)
        return missing
    end
    A = pi * (inner_diameter(e)/2)^2
    velocity = e.mass_flow / (WATER_DENSITY * A)
    return velocity
end

"""Compute water velocities for all pipe edges in the network.

Returns a dictionary keyed by `(src_label, dst_label)`.
"""
function water_velocities(nw::Network)::Dict{Tuple{String, String}, Float64}
    velocities = Dict{Tuple{String, String}, Float64}()
    for e in edges(nw)
        edge = nw[src(e), dst(e)]
        if edge isa InsulatedPipe
            velocities[(label_for(nw, src(e)), label_for(nw, dst(e)))] = water_velocity(edge)
        end
    end
    return velocities
end
