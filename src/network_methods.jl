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
Base.setindex!(nw::Network, v::SumpNode, label::String) = add_sump_node!(nw, v, label)
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
has_edge(nw::Network, src::String, dst::String) = Graphs.has_edge(nw.mg, index_for(nw, src), index_for(nw, dst))

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
        # check that sump nodes are internal (not leaves, not root)
        for sump_label in network.sump_labels
            if outdegree(network, sump_label) == 0
                error("Sump nodes must not be leaves (outdegree must be > 0)! Node $sump_label has outdegree 0.")
            end
            if indegree(network, sump_label) == 0
                error("Sump nodes must have exactly one incoming edge (indegree must be > 0)! Node $sump_label has indegree 0.")
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
    elseif v isa SumpNode
        # removing sump node
        filter!(x -> x != label, nw.sump_labels)
    end
    Graphs.rem_vertex!(nw.mg, index_for(nw, label))
    nw.neighbor_dicts.need_rebuild = true
end
rem_node!(nw::Network, idx::Integer) = rem_node!(nw, label_for(nw, idx))

"""Remove an edge from a `Network` by its source and destination string labels.

Both endpoint nodes remain in the network — only the connection between them is deleted.
`producer_label` and `load_labels` are not affected.
"""
function remove_edge!(nw::Network, src::String, dst::String)
    if !has_edge(nw, src, dst)
        error("Edge does not exist between $src and $dst")
    end
    Graphs.rem_edge!(nw.mg, index_for(nw, src), index_for(nw, dst))
    nw.neighbor_dicts.need_rebuild = true
end

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
"""Return the labels of all nodes that `label` has outgoing edges to."""
function outneighbors(nw::Network, label::String)
    check_and_update_neighbor_dicts!(nw)
    nw.neighbor_dicts.outneighbors[label]
end
"""Return the labels of all nodes that have edges pointing into `label`."""
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
position!(v::NodeType, pos::Tuple{Float64, Float64}) = (v.common.position = pos)

# vector operations for positions
Base.:+(a::Tuple{Float64, Float64}, b::Tuple{Float64, Float64}) = (a[1] + b[1], a[2] + b[2])
Base.:-(a::Tuple{Float64, Float64}, b::Tuple{Float64, Float64}) = (a[1] - b[1], a[2] - b[2])
Base.:*(a::Real, b::Tuple{Float64, Float64}) = (a * b[1], a * b[2])
Base.:*(a::Tuple{Float64, Float64}, b::Real) = (a[1] * b, a[2] * b)
Base.:/(a::Tuple{Float64, Float64}, b::Real) = (a[1] / b, a[2] / b)

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
function add_sump_node!(nw::Network, v::SumpNode, label::String)
    add_node!(nw, v, label)
    push!(nw.sump_labels, label)
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
        elseif nw[label] isa SumpNode
            # removing sump node from sump_labels
            filter!(x -> x != label, nw.sump_labels)
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

"""Check that a load function is non-negative and bounded over a temperature range.

Evaluates the function at a grid of sample points determined by the active flags and raises
an error if any value is **negative** or **not finite**.

# Arguments
- `fn`: load function. Signature depends on the flags (see [`LoadSpec`](@ref)).
- `params`: parameter vector passed to `fn`.
- `T_a_range`: ambient temperature sample points. Default: `-30.0:1.0:30.0`.
- `use_mass_flow`: whether `fn` takes mass flow as a third argument.
- `mass_flow_refs`: mass flow sample values used when `use_mass_flow = true`.
- `use_time`: whether `fn` takes a time-step index as a final integer argument.
- `time_refs`: integer time-step sample values used when `use_time = true`. Default: `1:10`.
"""
function validate_load_spec(fn::Function, params::Vector{Float64};
                             T_a_range=-30.0:1.0:30.0,
                             use_mass_flow::Bool=false,
                             mass_flow_refs::Vector{Float64}=[0.5, 1.0, 5.0],
                             use_time::Bool=false,
                             time_refs::AbstractVector{Int}=collect(1:10))
    if use_mass_flow && use_time
        for T_a in T_a_range, m in mass_flow_refs, t in time_refs
            v = fn(params, Float64(T_a), m, t)
            ismissing(v) && continue
            v < 0.0    && error("Load function returns negative power ($v kW) at T_ambient = $T_a °C, mass_flow = $m kg/s, time = $t. Power demand must be ≥ 0.")
            !isfinite(v) && error("Load function returns non-finite power ($v kW) at T_ambient = $T_a °C, mass_flow = $m kg/s, time = $t.")
        end
    elseif use_mass_flow
        for T_a in T_a_range, m in mass_flow_refs
            v = fn(params, Float64(T_a), m)
            ismissing(v) && continue
            v < 0.0    && error("Load function returns negative power ($v kW) at T_ambient = $T_a °C, mass_flow = $m kg/s. Power demand must be ≥ 0.")
            !isfinite(v) && error("Load function returns non-finite power ($v kW) at T_ambient = $T_a °C, mass_flow = $m kg/s.")
        end
    elseif use_time
        for T_a in T_a_range, t in time_refs
            v = fn(params, Float64(T_a), t)
            ismissing(v) && continue
            v < 0.0    && error("Load function returns negative power ($v kW) at T_ambient = $T_a °C, time = $t. Power demand must be ≥ 0.")
            !isfinite(v) && error("Load function returns non-finite power ($v kW) at T_ambient = $T_a °C, time = $t.")
        end
    else
        for T_a in T_a_range
            v = fn(params, Float64(T_a))
            ismissing(v) && continue
            v < 0.0    && error("Load function returns negative power ($v kW) at T_ambient = $T_a °C. Power demand must be ≥ 0.")
            !isfinite(v) && error("Load function returns non-finite power ($v kW) at T_ambient = $T_a °C.")
        end
    end
end

"""Set the power demand function and parameters for a single load node.

Validates the function over `T_a_range` before storing. The function must return power in **kW**
and be non-negative and finite over the validation range. The expected signature depends on the
flags — see [`LoadSpec`](@ref) and [`validate_load_spec`](@ref).

See also: [`set_load_params!`](@ref), [`validate_load_spec`](@ref).
"""
function set_load_fn!(nw::Network, label::String, fn::Function, params::AbstractVector{<:Real};
                      T_a_range=-30.0:1.0:30.0, use_mass_flow::Bool=false, use_time::Bool=false)
    has_label(nw, label) || error("Node with label $label does not exist in the network.")
    node = nw[label]
    node isa LoadNode || error("Node with label $label is not a load node.")
    p = Vector{Float64}(params)  # always copy so each node owns its params
    try
        validate_load_spec(fn, p; T_a_range=T_a_range, use_mass_flow=use_mass_flow, use_time=use_time)
    catch e
        error("$label power validation failed: $(e.msg)")
    end
    node.load = LoadSpec(fn, p, use_mass_flow, use_time)
    nw[label] = node
end

"""Set the same power demand function on every load node in the network, with the same parameters.

See also: [`set_load_fn!`](@ref).
"""
function set_load_fn!(nw::Network, fn::Function, params::AbstractVector{<:Real};
                      T_a_range=-30.0:1.0:30.0, use_mass_flow::Bool=false, use_time::Bool=false)
    for label in nw.load_labels
        set_load_fn!(nw, label, fn, params; T_a_range=T_a_range, use_mass_flow=use_mass_flow, use_time=use_time)
    end
end

"""Set the same power demand function on every load node, with per-load parameters from a dictionary.

`params_dict` maps each load label to its own parameter vector. All entries are validated before any
node is updated — if any validation fails, no nodes are changed.

See also: [`set_load_fn!`](@ref).
"""
function set_load_fn!(nw::Network, fn::Function, params_dict::Dict{String, <:AbstractVector{<:Real}};
                      T_a_range=-30.0:1.0:30.0, use_mass_flow::Bool=false, use_time::Bool=false)
    # validate all first so we don't partially update the network on error
    for (label, params) in params_dict
        has_label(nw, label) || error("Node with label $label does not exist in the network.")
        nw[label] isa LoadNode || error("Node with label $label is not a load node.")
        validate_load_spec(fn, convert(Vector{Float64}, params); T_a_range=T_a_range, use_mass_flow=use_mass_flow, use_time=use_time)
    end
    for (label, params) in params_dict
        node = nw[label]
        node.load = LoadSpec(fn, convert(Vector{Float64}, params), use_mass_flow, use_time)
        nw[label] = node
    end
end


function set_load_spec!(nw::Network, label::String, spec::LoadSpec)
    has_label(nw, label) || error("Node with label $label does not exist in the network.")
    node = nw[label]
    node isa LoadNode || error("Node with label $label is not a load node.")
    validate_load_spec(spec.fn, spec.params; T_a_range=-30.0:1.0:30.0, use_mass_flow=spec.use_mass_flow, use_time=spec.use_time)
    node.load = spec
    nw[label] = node
end

"""Set the relative mass-flow coefficient `m_rel` for a single load node.

`m_rel` is a dimensionless fraction relative to the total producer mass flow.
It is used by [`steady_state_hydrodynamics!`](@ref) to distribute absolute flow
among all loads; a value of 1.0 means *average* load (when all loads sum to
`n_loads`).  In practice, normalise so that `mean(m_rel) = 1.0` across all loads.

Two forms are accepted:

- `m_rel::Float64` — constant fraction, used identically at every time step.
- `m_rel::Vector{Float64}` — time-varying fraction; element `i` is used at simulation
  step `i`. The vector length must equal the number of simulation time steps `N`.
  All loads in a network must use the same form (constant or time-varying).

See also: [`set_load_fn!`](@ref), [`set_load_params!`](@ref), [`steady_state_hydrodynamics!`](@ref).
"""
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

function set_load_m_rel!(nw::Network, label::String, m_rel::Vector{Float64})
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

function set_load_m_rel!(nw::Network, m_rel::Dict{String, Float64})
    if Set(keys(m_rel)) != Set(nw.load_labels)
        throw(ArgumentError("Keys of m_rel dictionary must match load_labels exactly. Expected keys: $(nw.load_labels), got: $(keys(m_rel))"))
    end
    
    for label in nw.load_labels
        set_load_m_rel!(nw, label, m_rel[label])
    end
end

# ------------------------------------------------ #
# EDGE TYPE SPECIFIC METHODS
# ------------------------------------------------ #

"""Return physical pipe length in meters."""
pipe_length(e::InsulatedPipe) = e.physical_params.length

"""Same as [`pipe_length`](@ref) — returns the pipe length in meters."""
Base.length(e::InsulatedPipe) = pipe_length(e)

"""Return pipe inner diameter in meters."""
inner_diameter(e::InsulatedPipe) = e.physical_params.inner_diameter

"""
    volume(e::InsulatedPipe) -> Float64
    volume(e::ZeroPipe) -> Float64

Return the internal water volume of pipe `e` in cubic meters.

For an `InsulatedPipe` this is computed as `π/4 · L · d²`, where `L` is the pipe length
and `d` is the inner diameter. For a `ZeroPipe` the volume is always `0.0`.

# Examples
```julia
params = PipeParams(length=100.0, inner_diameter=0.05)
pipe = InsulatedPipe(info="supply", physical_params=params)
volume(pipe)   # ≈ 0.196 m³  (π/4 × 100 × 0.05²)
```
"""
volume(e::InsulatedPipe) = π/4 * pipe_length(e) * inner_diameter(e)^2
volume(e::ZeroPipe) = 0.0

"""Return thermal resistance (supply direction) in m·K/W."""
heat_resistance_forward(e::InsulatedPipe) = e.physical_params.heat_resistance_forward

"""Return thermal resistance (return direction) in m·K/W."""
heat_resistance_backward(e::InsulatedPipe) = e.physical_params.heat_resistance_backward

heat_resistance_forward!(e::InsulatedPipe, R::Float64) = (e.physical_params.heat_resistance_forward = R)
heat_resistance_backward!(e::InsulatedPipe, R::Float64) = (e.physical_params.heat_resistance_backward = R)


"""Return mass flow in kg/s."""
mass_flow(e::InsulatedPipe) = e.mass_flow

"""Return the steady-state mass flow [kg/s] stored on a node (set by `steady_state_hydrodynamics!`)."""
mass_flow(n::NodeType) = n.common.mass_flow

"""Return relative mass flow coefficient (`m_rel`) for a pipe or node.

    m_rel(e)           -> Union{Missing, Float64, Vector{Float64}}
    m_rel(e, step::Int) -> Float64

The single-argument form returns the raw field value — `missing` before assignment,
`Float64` after a per-step [`set_relative_mass_flows!`](@ref) call, or
`Vector{Float64}` after the vectorised no-argument overload.

The two-argument form always returns a `Float64` for the given `step`:
- `Float64` field → returns it directly (step is ignored; constant coefficient).
- `Vector{Float64}` field → returns `pipe.m_rel[step]`.

For a `LoadNode` the behaviour is the same: `m_rel(node)` returns
`Union{Missing, Float64, Vector{Float64}}`, and `m_rel(node, step)` returns
the scalar for that step.

See also: [`set_m_rel!`](@ref), [`set_load_m_rel!`](@ref),
[`set_relative_mass_flows!`](@ref).
"""
m_rel(e::InsulatedPipe) = e.m_rel
m_rel(n::NodeType) = n.m_rel

# Internal helper: extract the scalar m_rel value for a pipe at the given step.
_pipe_m_rel_at(m::Float64,          ::Int) = m
_pipe_m_rel_at(m::Vector{Float64}, step::Int) = length(m) == 1 ? m[1] : m[step]

m_rel(e::InsulatedPipe, step::Int) = _pipe_m_rel_at(e.m_rel, step)

"""Set the `m_rel` value on a pipe edge in-place.

    set_m_rel!(e::Union{InsulatedPipe, ZeroPipe}, value::Float64)
    set_m_rel!(e::Union{InsulatedPipe, ZeroPipe}, value::Vector{Float64})

Writes `value` to the `m_rel` field of the pipe. Both a constant scalar and a
pre-computed time-varying vector are accepted.

This is the low-level setter used internally by [`set_relative_mass_flows!`](@ref)
and is exposed for users who want to manipulate pipe split coefficients directly
(e.g. in custom hydraulic solvers).

After storing a `Vector{Float64}`, use [`m_rel(pipe, step)`](@ref) to read the
value for a specific time step.

See also: [`m_rel`](@ref), [`set_load_m_rel!`](@ref), [`set_relative_mass_flows!`](@ref).
"""
function set_m_rel!(e::Union{InsulatedPipe, ZeroPipe}, value::Union{Float64, Vector{Float64}})
    e.m_rel = value
end
function set_m_rel!(n::NodeType, value::Union{Float64, Vector{Float64}})
    n.m_rel = value
end

info(n::NodeType) = n.common.info
info(e::Union{InsulatedPipe, ZeroPipe}) = e.info
info!(n::NodeType, new_info::String) = (n.common.info = new_info)
info!(e::Union{InsulatedPipe, ZeroPipe}, new_info::String) = (e.info = new_info)

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

pipe_length(::ZeroPipe) = 0.0
Base.length(::ZeroPipe) = 0.0
inner_diameter(::ZeroPipe) = 0.0
heat_resistance_forward(::ZeroPipe) = 0.0
heat_resistance_backward(::ZeroPipe) = 0.0
mass_flow(e::ZeroPipe) = e.mass_flow
m_rel(e::ZeroPipe) = e.m_rel
m_rel(e::ZeroPipe, step::Int) = _pipe_m_rel_at(e.m_rel, step)
water_velocity(::ZeroPipe) = missing

"""Compute water velocities for all pipe edges in the network.

Returns a `Dict{Tuple{String,String}, Float64}` keyed by `(src_label, dst_label)`, with velocity in m/s.
Requires mass flows to be set first (e.g. via `steady_state_hydrodynamics!`).
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
