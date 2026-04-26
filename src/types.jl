# ------------------------------------------------ #
# TYPE DEFINITIONS
# ------------------------------------------------ #


"""Abstract supertype for all node types in a district heating network.

Concrete node types describe the *role* of a vertex in the directed network:

- `ProducerNode`: heat source (root)
- `JunctionNode`: branching/merging point
- `SumpNode`: measurement point (behaves like a junction, but recorded in `SimulationResults`)
- `LoadNode`: heat consumer (typically a leaf)
- `EmptyNode`: placeholder used during construction

All concrete node types store a `common::NodeCommon` field with metadata such as `info`, `position`, and (optionally) the current
steady-state `mass_flow`.
"""
abstract type NodeType end

"""Abstract supertype for all edge types in a district heating network.

Edges represent physical connections between nodes. In this package the primary edge type is `InsulatedPipe`.
`EmptyEdge` is used as a placeholder (e.g., when building a network from a bare topology).
"""
abstract type EdgeType end


"""Common data for all node types in the DH network.

`NodeCommon` contains fields that are useful across producers, junctions, and loads.

```julia
mutable struct NodeCommon
    info::String
    position::Union{Missing, Tuple{Float64, Float64}}   # (x, y) coordinates
    mass_flow::Union{Missing, Float64}                  # [kg/s]
end
```

# Fields
- `info::String`: human-readable label shown in printing/plots.
- `position::Union{Missing, Tuple{Float64, Float64}}`: optional `(x,y)` coordinates (used by visualization).
- `mass_flow::Union{Missing, Float64}`: steady-state mass flow through the node in kg/s.

# Notes
- Many fields use `missing` to represent “not initialized / not computed yet”.
- During simulation, mass flows are usually filled by `steady_state_hydrodynamics!`.

# Constructors
- `NodeCommon(info::String)`
- `NodeCommon(info::String, position::Tuple{Float64, Float64})`
"""
mutable struct NodeCommon
    info::String
    position::Union{Missing, Tuple{Float64, Float64}}   # (x, y) coordinates
    mass_flow::Union{Missing, Float64}                  # [kg/s]
end


"""DH network node representing a junction.

Junctions are internal vertices that connect multiple pipes but do not directly produce or consume heat.
They are where flow splits (supply direction) or merges (return direction).

```julia
struct JunctionNode <: NodeType
    common::NodeCommon
end
```

# Constructors
- `JunctionNode(info::String)`
- `JunctionNode(info::String, position::Tuple{Float64, Float64})`
- `JunctionNode(position::Tuple{Float64, Float64})`
- `JunctionNode()`
"""
struct JunctionNode <: NodeType
    common::NodeCommon
end

"""DH network node representing a sump — a measurement point in the network.

A `SumpNode` behaves identically to a [`JunctionNode`](@ref) during simulation: it routes
flow between pipes without consuming or producing heat. The difference is that sumps are
*tracked* — supply temperature, return temperature, and total mass flow at every sump are
recorded in [`SimulationResults`](@ref) at each time step.

Use sumps wherever you want to observe network conditions at an intermediate node without
disturbing the topology.

```julia
struct SumpNode <: NodeType
    common::NodeCommon
end
```

# Constructors
- `SumpNode(info::String)`
- `SumpNode(info::String, position::Tuple{Float64, Float64})`
- `SumpNode(position::Tuple{Float64, Float64})`
- `SumpNode()`

See also: [`JunctionNode`](@ref), [`SimulationResults`](@ref).
"""
struct SumpNode <: NodeType
    common::NodeCommon
end

"""Built-in power demand function: polynomial in ambient temperature.

Evaluates ``P(T_a) = \\sum_{i} params_i \\cdot T_a^{i-1}`` and returns power in **kW**.
Any number of coefficients is supported:

- `params = [p₀, p₁, p₂]` — quadratic
- `params = [p₀, p₁]` — linear
- `params = [p₀]` — constant

The function returns zero for ambient temperatures above 30 °C to prevent unrealistic behavior in the summer.

See also: [`LoadSpec`](@ref).
"""
polynomial_load(params::Vector{Float64}, T_a::Float64) = T_a < 30 ? sum(params[i] * T_a^(i-1) for i in eachindex(params)) : 0

"""Built-in power demand function: two-part (hockey-stick) linear model.

Models heating demand as a piecewise-linear function of ambient temperature:

```
P(T_a) = a + b·(T_b − T_a)   if T_a ≤ T_b   (heating regime)
P(T_a) = a                    if T_a >  T_b   (base-load regime)
```

The result is clamped to zero from below so it never goes negative.

# Parameters (`params = [a, b, T_b]`)
- `a`  — base load [kW]: the constant demand present at all temperatures.
- `b`  — heating slope [kW/°C]: additional load per degree below `T_b`.
- `T_b` — balance-point temperature [°C]: above this the load is flat at `a`.

# Example
```julia
# 50 kW base load, 5 kW per °C below 15 °C
spec = LoadSpec(hockey_load, [50.0, 5.0, 15.0])
```

See also: [`LoadSpec`](@ref), [`polynomial_load`](@ref).
"""
function hockey_load(params::Vector{Float64}, T_a::Float64)
    a, b, T_b = params[1], params[2], params[3]
    @assert b >= 0 "Heating slope b must be non-negative"
    @assert a >= 0 "Base load a must be non-negative"
    return max(0.0, T_a <= T_b ? a + b * (T_b - T_a) : a)
end

"""Built-in power demand function: general hockey-stick model depending on ambient
temperature **and** mass flow through the load.

The load parameters `a` and `b` of the hockey-stick model are expressed as linear functions
of the mass flow ṁ [kg/s] at the load:

```
a(ṁ)      = α_a + β_a · ṁ
b(ṁ)      = α_b + β_b · ṁ
P(T_a, ṁ) = a(ṁ) + b(ṁ) · max(0, T_b − T_a)   [kW]
```

# Parameters (`params = [α_a, β_a, α_b, β_b, T_b]`)
- `α_a` — base-load intercept [kW].
- `β_a` — base-load slope [kW/(kg/s)].
- `α_b` — heating-slope intercept [kW/°C].
- `β_b` — heating-slope slope [kW/(°C·kg/s)].
- `T_b` — balance-point temperature [°C].

See also: [`hockey_load`](@ref), [`LoadSpec`](@ref).
"""
function general_hockey_load(params::Vector{Float64}, T_a::Float64, mass_flow::Float64)
    α_a, β_a, α_b, β_b, T_b = params
    a = α_a + β_a * mass_flow
    b = α_b + β_b * mass_flow
    return max(0.0, a + b * max(0.0, T_b - T_a))
end

"""Load specification pairing a power demand function with its parameters.

The load function can depend on ambient temperature alone, or on both ambient temperature
and the load's instantaneous mass flow (e.g. [`general_hockey_load`](@ref)).

```julia
mutable struct LoadSpec
    fn::Function
    params::Vector{Float64}
    use_mass_flow::Bool
end
```

# Fields
- `fn::Function`: demand function. When `use_mass_flow = false`, called as
  `fn(params, T_a) -> Float64` [kW]. When `use_mass_flow = true`, called as
  `fn(params, T_a, mass_flow) -> Float64` [kW], where `mass_flow` is the
  steady-state mass flow [kg/s] at the load node.
- `params::Vector{Float64}`: current parameters passed to `fn`.
- `use_mass_flow::Bool`: whether the function takes mass flow as a third argument.

See also: [`polynomial_load`](@ref), [`hockey_load`](@ref), [`general_hockey_load`](@ref),
[`set_load_fn!`](@ref), [`set_load_params!`](@ref), [`validate_load_spec`](@ref).

# Constructors
- `LoadSpec()`: default `polynomial_load` with `DEFAULT_LOAD_PARAMS`.
- `LoadSpec(params::Vector{<:Real})`: `polynomial_load` with custom parameter vector.
- `LoadSpec(p₀, p₁, ...)`: `polynomial_load` with parameters given as individual numbers.
- `LoadSpec(fn, params::Vector{<:Real})`: explicit function + parameter vector (converts element type).
- `LoadSpec(fn, p₀, p₁, ...)`: explicit function with parameters as individual numbers.
- `LoadSpec(hockey_load; a, b, T_b)`: keyword form for `hockey_load` — base load, slope, balance-point temperature.
"""
mutable struct LoadSpec
    fn::Function
    params::Vector{Float64}
    use_mass_flow::Bool
end

# Load Node: heat consumer
"""DH network node representing a load (consumer).

The `load` field holds a [`LoadSpec`](@ref) pairing a demand function with its parameters.
The function is called as `load.fn(load.params, T_a)` and must return power demand in **kW**
for a given ambient temperature `T_a` in °C.

`m_rel` controls how the total flow is divided between branches. A load with `m_rel=2.0` gets twice as much flow as
one with `m_rel=1.0`. Set this on every load node before running a simulation — the solver reads it from the leaves
and propagates the split ratios upstream.

`m_rel` can be either a **constant** `Float64` (same split every time step) or a **time-varying** `Vector{Float64}`
(one value per simulation time step). All loads in the network must use the same mode uniformly: mixing constant and
time-varying loads is not allowed and is detected at the start of [`run_simulation`](@ref).

When a `Vector{Float64}` is stored, its length must equal the number of time steps `N` passed to [`run_simulation`](@ref).

```julia
mutable struct LoadNode <: NodeType
    common::NodeCommon
    load::Union{Missing, LoadSpec}                      # Power demand specification
    m_rel::Union{Missing, Float64, Vector{Float64}}     # Relative mass flow coefficient (constant or time-varying)
end
```

# Constructors
- `LoadNode(; load=LoadSpec(polynomial_load, [...]))`
- `LoadNode(info::String; load=...)`
- `LoadNode(info::String, m_rel::Float64; load=...)`: constant `m_rel`, no position.
- `LoadNode(info::String, m_rel::Vector{Float64}; load=...)`: time-varying `m_rel`, no position.
- `LoadNode(position::Tuple{Float64, Float64}; load=...)`
- `LoadNode(info::String, position::Tuple{Float64, Float64}; load=...)`
- `LoadNode(info::String, position::Tuple{Float64, Float64}, m_rel::Float64; load=...)`: constant `m_rel` at construction.
- `LoadNode(info::String, position::Tuple{Float64, Float64}, m_rel::Vector{Float64}; load=...)`: time-varying `m_rel` at construction.
- `LoadNode(info::String, position::Tuple{Float64, Float64}, load::LoadSpec)`: provide a custom `LoadSpec` as the third positional argument.

See also: [`set_load_m_rel!`](@ref), [`run_simulation`](@ref).
"""
mutable struct LoadNode <: NodeType
    common::NodeCommon
    load::Union{Missing, LoadSpec}                      # Power demand specification
    m_rel::Union{Missing, Float64, Vector{Float64}}     # Relative mass flow coefficient (constant or time-varying)
end

# Producer Node: heat producer
"""DH network node representing a producer (heat source).

There may be only one producer in the network, and it is identified by its label (producer_label field in Network).

In simulations, the producer’s mass flow and supply temperature are usually provided by a control policy
(see `ProducerOutput` and `run_simulation`).

```julia
struct ProducerNode <: NodeType
    common::NodeCommon
end
```

# Constructors
- `ProducerNode(info::String)`
- `ProducerNode(info::String, position::Tuple{Float64, Float64})`
- `ProducerNode(position::Tuple{Float64, Float64})`
- `ProducerNode()`
"""
struct ProducerNode <: NodeType
    common::NodeCommon
end

# empty Node: just for initialization purposes
"""Placeholder node type.

`EmptyNode` is used when constructing a `Network` from a bare topology (e.g. `Network(graph::DiGraph)`).
Replace placeholders with real node types (`ProducerNode`, `JunctionNode`, `LoadNode`) before running simulations.
"""
struct EmptyNode <: NodeType end

# ------------------------------------------------ #
# PIPE MODEL
# ------------------------------------------------ #
# we will model pipe using "plug" method

"""A single water plug used by the plug-method pipe model.

A `Plug` represents a mass of water inside a pipe that is assumed to have uniform temperature.
Plugs are advected through pipes during time stepping and may be split/merged.

```julia
mutable struct Plug
    T::Float64  # Temperature at the plug [°C]
    m::Float64  # mass of the plug [kg]
    k::Float64  # fractional step when the plug's midpoint entered the current pipe (0.0 = initial fill)
end
```

`k` is a fractional (sub-step) entry time in units of simulation steps. When a plug exits a
pipe, `k` is updated to the fractional step at which the plug's mass midpoint crossed the pipe
outlet — this becomes the entry time for the next pipe segment. Transit time through any pipe is
`τ = (k_exit - k_entry) · Δt`, giving sub-step resolution without integer rounding.
"""
mutable struct Plug
    T::Float64  # Temperature at the plug [°C]
    m::Float64  # mass of the plug [kg]
    k::Float64  # fractional step when the plug's midpoint entered the current pipe (0.0 = initial fill)
end

"""    Plug(T, m)
Convenience constructor — creates a plug with `k = 0.0` (initial fill, not yet in a pipe).
"""
Plug(T::Float64, m::Float64) = Plug(T, m, 0.0)

"""Physical parameters of a pipe (constant during simulation).

These parameters describe geometry and heat-loss characteristics.

```julia
mutable struct PipeParams                       # unchanging physical parameters of the pipe
    length::Float64                     # Length of the pipe [m]
    inner_diameter::Float64             # Inner diameter [m]
    heat_resistance_forward::Float64    # Thermal resistance [m*K/W]
    heat_resistance_backward::Float64   # Thermal resistance [m*K/W]
end
```

# Fields
- `length` [m]
- `inner_diameter` [m]
- `heat_resistance_forward` [m·K/W]: thermal resistance for the supply direction
- `heat_resistance_backward` [m·K/W]: thermal resistance for the return direction

# Constructors
- `PipeParams(length, inner_diameter, heat_resistance_forward, heat_resistance_backward)`: positional, all values explicit.
- `PipeParams(; length, inner_diameter, heat_resistance_forward=3.0, heat_resistance_backward=4.0)`: keyword form with sensible defaults for the resistance values.
"""
mutable struct PipeParams                       # unchanging physical parameters of the pipe
    length::Float64                     # Length of the pipe [m]
    inner_diameter::Float64             # Inner diameter [m]
    heat_resistance_forward::Float64    # Thermal resistance [m*K/W]
    heat_resistance_backward::Float64   # Thermal resistance [m*K/W]
end

# InsulatedPipe: represents a pipe in the network
"""DH network edge representing an insulated pipe.

`InsulatedPipe` transports water between nodes and stores both physical parameters and the current hydraulic/thermal state.

```julia
mutable struct InsulatedPipe <: EdgeType
    info::String
    physical_params::PipeParams
    mass_flow::Union{Missing, Float64}                  # Mass flow in [kg/s]
    m_rel::Union{Missing, Float64, Vector{Float64}}     # Relative mass flow coefficient (scalar or time-varying)
    plugs_f::Vector{Plug}               # Queue of plugs in the pipe (forward direction)
    plugs_b::Vector{Plug}               # Queue of plugs in the pipe (backward direction)
end
```

# Fields
- `physical_params`: geometry and heat-loss parameters (see `PipeParams`).
- `mass_flow`: mass flow in kg/s. `missing` until computed by `steady_state_hydrodynamics!`.
- `m_rel`: relative flow coefficient used for splitting at junctions. `missing` until set.
  Written as a `Float64` scalar by the per-step [`set_relative_mass_flows!`](@ref) overload
  or as a `Vector{Float64}` (one entry per simulation step) by the no-argument (vectorised)
  overload. Use [`m_rel(pipe, step)`](@ref) to read the value for a specific time step
  regardless of which form is stored.
- `plugs_f`: plug queue for the forward (supply) direction.
- `plugs_b`: plug queue for the backward (return) direction.

# Constructors
- `InsulatedPipe(; info="pipe", length=100.0, inner_diameter=0.1, heat_resistance_forward=3.0, heat_resistance_backward=4.0, mass_flow=missing, m_rel=missing)`: all-keyword constructor using individual pipe dimensions; physical params have sensible defaults.
- `InsulatedPipe(; info="pipe", params::PipeParams, mass_flow=missing, m_rel=missing)`: all-keyword constructor accepting a pre-constructed `PipeParams`; individual dimension keywords are ignored when `params` is supplied.
- `InsulatedPipe(info::String, params::PipeParams; mass_flow=missing, m_rel=missing)`: positional-style build from a pre-constructed `PipeParams`.
- `InsulatedPipe(params::PipeParams; mass_flow=missing, m_rel=missing)`: same, with default `info="pipe"`.
- `InsulatedPipe(length::Real)`: shorthand that sets only the pipe length; all other params take their defaults.
"""
mutable struct InsulatedPipe <: EdgeType
    info::String
    physical_params::PipeParams
    mass_flow::Union{Missing, Float64}                  # Mass flow in [kg/s]
    m_rel::Union{Missing, Float64, Vector{Float64}}     # Relative mass flow coefficient (scalar or time-varying)
    plugs_f::Vector{Plug}               # Queue of plugs in the pipe (forward direction)
    plugs_b::Vector{Plug}               # Queue of plugs in the pipe (backward direction)
end

"""Placeholder edge type.

`EmptyEdge` is used when constructing a `Network` from a topology without pipe parameters.
Replace it with `InsulatedPipe` before running simulations.
"""
struct EmptyEdge <: EdgeType end

# ------------------------------------------------- #
# DH NETWORK TYPE
# ------------------------------------------------- #

"""Cached neighbor lookup tables stored inside a `Network`.

Instead of scanning the graph every time `outneighbors` or `inneighbors` is called, the results
are pre-computed once and stored here. Since the network topology is fixed during a simulation,
this cache is built once before the first time step and reused throughout.

The `need_rebuild` flag is set to `true` whenever the topology changes (e.g. after adding or
removing nodes/edges), so the cache is automatically refreshed before the next lookup.

# Constructor
- `NeighborDicts()`: creates an empty instance with `need_rebuild = true`.
"""
mutable struct NeighborDicts
    outneighbors::Dict{String, Vector{String}}       # mapping nodes to outneighbors for efficient access in simulation
    inneighbors::Dict{String, Vector{String}}        # mapping nodes to inneighbors for efficient access in simulation
    need_rebuild::Bool                               # flag to indicate if neighbor dicts need to be rebuilt before simulation
end
NeighborDicts() = NeighborDicts(Dict{String, Vector{String}}(), Dict{String, Vector{String}}(), true)


"""Network type representing a district heating network.
```julia
mutable struct Network{T<:Integer} <: AbstractGraph{T}
    mg::MetaGraph                               # MetaGraph from MetaGraphs.jl, it contains the network
                                                    structure and node and edge data
    producer_label::Union{Nothing, String}      # Label of the producer node
    load_labels::Set{String}                    # Labels of the consumer nodes
    sump_labels::Set{String}                    # Labels of the sump nodes
    neighbor_dicts::NeighborDicts               # Mappings to inneghbors and outneighbors for efficient access during simulation
end
```

# Fields
- `mg`: a MetaGraphsNext `MetaGraph` that contains the directed topology and stores node/edge data.
- `producer_label`: label of the single producer node (or `nothing` if not yet set).
- `load_labels`: a set of labels for load nodes.
- `sump_labels`: a set of labels for sump nodes.
- `neighbor_dicts`: cached neighbor lists used to reduce allocations during simulation.

# Constructors
- `Network()`: creates an empty network with no nodes or edges.
- `Network(g::SimpleDiGraph)`: wraps an existing `Graphs.jl` directed graph. All nodes and edges are
  initialized as `EmptyNode` / `EmptyEdge` — use `name_nodes!`, `identify_producer_and_loads!`, etc.
  to populate them.

"""
mutable struct Network{T<:Integer} <: AbstractGraph{T}
    mg::MetaGraph                               # MetaGraph from MetaGraphs.jl, it contains the network structure and node and edge data
    producer_label::Union{Nothing, String}      # Label of the producer node
    load_labels::Set{String}                    # Labels of the consumer nodes
    sump_labels::Set{String}                    # Labels of the sump nodes
    neighbor_dicts::NeighborDicts               # Mappings to inneghbors and outneighbors for efficient access during simulation

    function Network(mg::MetaGraph, producer_label, load_labels, sump_labels, neighbor_dicts)
        network = new{eltype(vertices(mg))}(mg, producer_label, load_labels, sump_labels, neighbor_dicts)
        check_and_update_neighbor_dicts!(network) # upon creation, construct the neighbor dicts
        return network
    end
end

# ------------------------------------------------- #
# CONSTRUCTORS
# ------------------------------------------------- #

# default constructor for empty network, more about construction in network_creation.jl
function Network()
    mg = MetaGraph(
        DiGraph();  # underlying graph structure
        label_type=String,
        vertex_data_type=NodeType,
        edge_data_type=EdgeType
    )
    return Network(mg, nothing, Set{String}(), Set{String}(), NeighborDicts())
end

# ------------------------------------------------- #
# PIPE CONSTRUCTORS
# ------------------------------------------------- #

PipeParams(;length::Real, inner_diameter::Real, heat_resistance_forward::Real=3.0, heat_resistance_backward::Real=4.0) = PipeParams(float(length), float(inner_diameter), float(heat_resistance_forward), float(heat_resistance_backward))

function InsulatedPipe(; info::String="pipe",
                         params::Union{PipeParams, Nothing}=nothing,
                         length::Real=100.0,
                         inner_diameter::Real=0.1,
                         heat_resistance_forward::Real=3.0,
                         heat_resistance_backward::Real=4.0,
                         mass_flow=missing,
                         m_rel=missing)
    if params === nothing
        params = PipeParams(float(length), float(inner_diameter),
                            float(heat_resistance_forward), float(heat_resistance_backward))
    end
    return InsulatedPipe(info, params, mass_flow, m_rel, Vector{Plug}(), Vector{Plug}())
end

function InsulatedPipe(info::String, params::PipeParams; mass_flow=missing, m_rel=missing)
    return InsulatedPipe(info, params, mass_flow, m_rel, Vector{Plug}(), Vector{Plug}())
end

InsulatedPipe(params::PipeParams; mass_flow=missing, m_rel=missing) =
    InsulatedPipe("pipe", params; mass_flow=mass_flow, m_rel=m_rel)

InsulatedPipe(length::Real) = InsulatedPipe(; length=float(length))

"""DH network edge representing an instant, lossless connection.

`ZeroPipe` models a direct connection between two nodes with no transit delay,
no heat loss, and no pressure drop. Use it wherever a real pipe is not needed
(e.g. the zero-delay branch of an aggregated junction pair).

Unlike `InsulatedPipe`, `ZeroPipe` carries no `physical_params`. Its length,
inner diameter, and thermal resistances are identically zero by definition.
The `plugs_f`/`plugs_b` queues act as within-timestep buffers — they are always
empty at the start and end of each time step.

# Constructors
- `ZeroPipe()`: default info `"zero pipe"`.
- `ZeroPipe(info::String)`: custom label.
- `ZeroPipe(info::String; mass_flow=missing, m_rel=missing)`: pre-set hydraulic fields.
- `ZeroPipe(; info="zero pipe", mass_flow=missing, m_rel=missing)`: keyword form.

`m_rel` follows the same scalar / time-varying convention as [`InsulatedPipe`](@ref): a
`Float64` after a per-step [`set_relative_mass_flows!`](@ref) call, or a `Vector{Float64}`
after the vectorised no-argument overload. Use [`m_rel(pipe, step)`](@ref) to read the
value for a specific time step.

See also: `InsulatedPipe`.
"""
mutable struct ZeroPipe <: EdgeType
    info::String
    mass_flow::Union{Missing, Float64}
    m_rel::Union{Missing, Float64, Vector{Float64}}
    plugs_f::Vector{Plug}
    plugs_b::Vector{Plug}
end

ZeroPipe(info::String; mass_flow=missing, m_rel=missing) = ZeroPipe(info, mass_flow, m_rel, Vector{Plug}(), Vector{Plug}())
ZeroPipe(; info::String="zero pipe", mass_flow=missing, m_rel=missing) = ZeroPipe(info; mass_flow=mass_flow, m_rel=m_rel)


# ------------------------------------------------- #
# NODE CONSTRUCTORS
# ------------------------------------------------- #

# NodeCommon constructors
NodeCommon(info::String) = NodeCommon(info, missing, missing)
NodeCommon(info::String, position::Tuple{Float64, Float64}) = NodeCommon(info, position, missing)

# JUNCTION NODE CONSTRUCTORS
JunctionNode(info::String) = JunctionNode(NodeCommon(info))
JunctionNode(info::String, position::Tuple{Float64, Float64}) = JunctionNode(NodeCommon(info, position))
JunctionNode(position::Tuple{Float64, Float64}) = JunctionNode("junction", position)
JunctionNode() = JunctionNode("junction")

# SUMP NODE CONSTRUCTORS
SumpNode(info::String) = SumpNode(NodeCommon(info))
SumpNode(info::String, position::Tuple{Float64, Float64}) = SumpNode(NodeCommon(info, position))
SumpNode(position::Tuple{Float64, Float64}) = SumpNode("sump", position)
SumpNode() = SumpNode("sump")

# LOAD NODE CONSTRUCTORS
const DEFAULT_LOAD_PARAMS = [540.0, -36.0, 0.6]  # default quadratic polynomial coefficients (kW vs °C)
LoadSpec() = LoadSpec(polynomial_load, copy(DEFAULT_LOAD_PARAMS), false)
LoadSpec(params::Vector{<:Real})       = LoadSpec(polynomial_load, Vector{Float64}(params), false)
LoadSpec(params::Real...)              = LoadSpec(polynomial_load, collect(Float64, params), false)
LoadSpec(fn::Function, params::Vector{<:Real}) = LoadSpec(fn, Vector{Float64}(params), false)
LoadSpec(fn::Function, params::Real...) = LoadSpec(fn, collect(Float64, params), false)
LoadSpec(::typeof(hockey_load); a::Real, b::Real, T_b::Real) = LoadSpec(hockey_load, [Float64(a), Float64(b), Float64(T_b)], false)

LoadNode(info::String; load::LoadSpec=LoadSpec()) = LoadNode(NodeCommon(info), load, missing)
LoadNode(info::String, m_rel::Float64; load::LoadSpec=LoadSpec()) = LoadNode(NodeCommon(info), load, m_rel)
LoadNode(info::String, m_rel::Vector{Float64}; load::LoadSpec=LoadSpec()) = LoadNode(NodeCommon(info), load, m_rel)
LoadNode(info::String, position::Tuple{Float64, Float64}; load::LoadSpec=LoadSpec()) = LoadNode(NodeCommon(info, position), load, missing)
LoadNode(info::String, position::Tuple{Float64, Float64}, load::LoadSpec) = LoadNode(NodeCommon(info, position), load, missing)
LoadNode(info::String, position::Tuple{Float64, Float64}, m_rel::Float64; load::LoadSpec=LoadSpec()) = LoadNode(NodeCommon(info, position), load, m_rel)
LoadNode(info::String, position::Tuple{Float64, Float64}, m_rel::Vector{Float64}; load::LoadSpec=LoadSpec()) = LoadNode(NodeCommon(info, position), load, m_rel)
LoadNode(info::String, ::Missing; load::LoadSpec=LoadSpec()) = LoadNode(info; load=load)  # position missing, use keyword form
LoadNode(position::Tuple{Float64, Float64}; load::LoadSpec=LoadSpec()) = LoadNode("load", position; load=load)
LoadNode(; load::LoadSpec=LoadSpec()) = LoadNode("load"; load=load)

# PRODUCER NODE CONSTRUCTORS
ProducerNode(info::String) = ProducerNode(NodeCommon(info))
ProducerNode(info::String, position::Tuple{Float64, Float64}) = ProducerNode(NodeCommon(info, position))
ProducerNode(position::Tuple{Float64, Float64}) = ProducerNode("producer", position)
ProducerNode() = ProducerNode("producer")
