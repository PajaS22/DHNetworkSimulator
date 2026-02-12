# ------------------------------------------------ #
# TYPE DEFINITIONS
# ------------------------------------------------ #


"""Abstract supertype for all node types in a district heating network.

Concrete node types describe the *role* of a vertex in the directed network:

- `ProducerNode`: heat source (root)
- `JunctionNode`: branching/merging point
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
- During simulation, mass flows are usually filled by `steady_state_hydronynamics!`.

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

# Load Node: heat consumer
"""DH network node representing a load (consumer).

The `load` field defines a quadratic power-demand curve as a function of ambient temperature ``T_a``:

``P(T_a) = p_0 + p_1 T_a + p_2 T_a^2,``
where the power is in kW and the ambient temperature is in °C.

The `m_rel` field is a *relative mass-flow coefficient* used when splitting flows at a junction while solving steady state flow;
it is specified for leaf nodes and propagated to upstream edges.

```julia
mutable struct LoadNode <: NodeType
    common::NodeCommon
    load::Union{Missing, NTuple{3, Float64}}     # Heat load function in kW, P(Tₐ) = p₀ + p₁*Tₐ + p₂*Tₐ²
    m_rel::Union{Missing, Float64}               # Relative mass flow coefficient (for branching nodes)
end
```

# Constructors
- `LoadNode(info::String; load=DEFAULT_LOAD)`: Creates a LoadNode with the specified info string and an optional load function (default is a typical load curve).
- `LoadNode(info::String, position::Tuple{Float64, Float64}; load=DEFAULT_LOAD)`
- `LoadNode(info::String, position::Tuple{Float64, Float64}, load::NTuple{3, Float64})`
- `LoadNode(info::String, position::Tuple{Float64, Float64}, m_rel::Float64; load=DEFAULT_LOAD)`
- `LoadNode(position::Tuple{Float64, Float64}; load=DEFAULT_LOAD)`
- `LoadNode(; load=DEFAULT_LOAD)`
"""
mutable struct LoadNode <: NodeType
    common::NodeCommon
    load::Union{Missing, NTuple{3, Float64}}     # Heat load function in kW, P(Tₐ) = p₀ + p₁*Tₐ + p₂*Tₐ²
    m_rel::Union{Missing, Float64}               # Relative mass flow coefficient (for branching nodes)
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
end
```
"""
mutable struct Plug
    T::Float64  # Temperature at the plug [°C]
    m::Float64  # mass of the plug [kg]
end

"""Physical parameters of a pipe (constant during simulation).

These parameters describe geometry and heat-loss characteristics.

```julia
struct PipeParams                       # unchanging physical parameters of the pipe
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
- `PipeParams(length::Float64, inner_diameter::Float64, heat_resistance_forward::Float64, heat_resistance_backward::Float64)`
- `PipeParams(length::Float64, inner_diameter::Float64)`: uses default heat resistance values based on typical insulation properties.
"""
struct PipeParams                       # unchanging physical parameters of the pipe
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
    mass_flow::Union{Missing, Float64}  # Mass flow in [kg/s]
    m_rel::Union{Missing, Float64}      # Relative mass flow coefficient (for branching pipes)
    plugs_f::Vector{Plug}               # Queue of plugs in the pipe (forward direction)
    plugs_b::Vector{Plug}               # Queue of plugs in the pipe (backward direction)
end
```

# Fields
- `physical_params: geometry and heat-loss parameters.
- `mass_flow: mass flow in kg/s (typically computed by `steady_state_hydronynamics!`).
- `m_rel: relative flow coefficient used for splitting at junctions.
- `plugs_f: plug queue for the forward (supply) direction.
- `plugs_b: plug queue for the backward (return) direction.

# Constructors
- `InsulatedPipe(info::String; length::Float64, inner_diameter::Float64, heat_resistance_forward::Float64, heat_resistance_backward::Float64)`
- `InsulatedPipe(info::String, params::PipeParams)`
- `InsulatedPipe(params::PipeParams)`
- `InsulatedPipe(length::Real)`
"""
mutable struct InsulatedPipe <: EdgeType
    info::String
    physical_params::PipeParams
    mass_flow::Union{Missing, Float64}  # Mass flow in [kg/s]
    m_rel::Union{Missing, Float64}      # Relative mass flow coefficient (for branching pipes)
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

"""Mappings to inneghbors and outneighbors for efficient access during simulation, stored in the Network struct.

When calling functions `outneighbors(nw, label)` or `inneighbors(nw, label)`, 
there is no need to scan through the graph and collect neighbors, 
instead we can directly access the pre-computed neighbor lists in the dictionaries.
This significantly lowers the number of allocations during simulation, because there
is a lot of places where we need to access neighbors of a node.

The network is static during the simulation, so we can compute these neighbor lists once
before the simulation starts and then reuse them.

```julia
mutable struct NeighborDicts
    outneighbors::Dict{String, Vector{String}}       # mapping nodes to outneighbors for efficient access in simulation
    inneighbors::Dict{String, Vector{String}}        # mapping nodes to inneighbors for efficient access in simulation
    need_rebuild::Bool                               # flag to indicate if neighbor dicts need to be rebuilt before simulation
end
```

- `need_rebuild` flag is used to indicate when the neighbor dicts need to be updated
   (e.g., after adding/removing nodes or edges), so that we can avoid unnecessary 
   rebuilding during multiple modifications.

# Constructor
- `NeighborDicts()`: Creates an instance of NeighborDicts with empty dictionaries and the need_rebuild flag set to true.
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
    neighbor_dicts::NeighborDicts               # Mappings to inneghbors and outneighbors for efficient access during simulation
end
```

# Fields
- `mg`: a MetaGraphsNext `MetaGraph` that contains the directed topology and stores node/edge data.
- `producer_label`: label of the single producer node (or `nothing` if not yet set).
- `load_labels`: a set of labels for load nodes.
- `neighbor_dicts`: cached neighbor lists used to reduce allocations during simulation.

# Constructors
- `Network()`: Creates an empty DH network with no nodes or edges.
- `Network(g::DiGraph)`: Creates a DH network from an existing directed graph structure ('Graphs.jl').
                         Nodes and edges are initialized with EmptyNode and EmptyEdge data.

"""
mutable struct Network{T<:Integer} <: AbstractGraph{T}
    mg::MetaGraph                               # MetaGraph from MetaGraphs.jl, it contains the network structure and node and edge data
    producer_label::Union{Nothing, String}      # Label of the producer node
    load_labels::Set{String}                    # Labels of the consumer nodes
    neighbor_dicts::NeighborDicts               # Mappings to inneghbors and outneighbors for efficient access during simulation
    
    function Network(mg::MetaGraph, producer_label, load_labels, neighbor_dicts) 
        network = new{eltype(vertices(mg))}(mg, producer_label, load_labels, neighbor_dicts)
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
    return Network(mg, nothing, Set{String}(), NeighborDicts())
end

# ------------------------------------------------- #
# PIPE CONSTRUCTORS
# ------------------------------------------------- #

function InsulatedPipe(info::String;
                        length::Float64=100.0,
                        inner_diameter::Float64=0.1,
                        heat_resistance_forward::Float64=3.0,
                        heat_resistance_backward::Float64=4.0)
    physical_params = PipeParams(length, inner_diameter, heat_resistance_forward, heat_resistance_backward)
    return InsulatedPipe(info, physical_params)
end
function InsulatedPipe(info::String, params::PipeParams)::InsulatedPipe
    mass_flow = missing         # [kg/s]
    m_rel = missing             # Relative mass flow coefficient
    plugs_f = Vector{Plug}()    # Initialize empty queue of plugs (forward direction)
    plugs_b = Vector{Plug}()    # Initialize empty queue of plugs (backward direction)
    return InsulatedPipe(info, params, mass_flow, m_rel, plugs_f, plugs_b)
end
InsulatedPipe(params::PipeParams) = InsulatedPipe("pipe", params)
InsulatedPipe(length::Real) = InsulatedPipe("pipe"; length=float(length))

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

# LOAD NODE CONSTRUCTORS
const DEFAULT_LOAD = (540.0, -36.0, 0.6)
LoadNode(info::String; load=DEFAULT_LOAD) = LoadNode(NodeCommon(info), load, missing)
LoadNode(info::String, position::Tuple{Float64, Float64}; load=DEFAULT_LOAD) = LoadNode(NodeCommon(info, position), load, missing)
LoadNode(info::String, position::Tuple{Float64, Float64}, load::NTuple{3, Float64}) = LoadNode(NodeCommon(info, position), load, missing)
LoadNode(info::String, position::Tuple{Float64, Float64}, m_rel::Float64; load=DEFAULT_LOAD) = LoadNode(NodeCommon(info, position), load, m_rel)
LoadNode(position::Tuple{Float64, Float64}; load=DEFAULT_LOAD) = LoadNode("load", position; load=load)
LoadNode(; load=DEFAULT_LOAD) = LoadNode("load"; load=load)

# PRODUCER NODE CONSTRUCTORS
ProducerNode(info::String) = ProducerNode(NodeCommon(info))
ProducerNode(info::String, position::Tuple{Float64, Float64}) = ProducerNode(NodeCommon(info, position))
ProducerNode(position::Tuple{Float64, Float64}) = ProducerNode("producer", position)
ProducerNode() = ProducerNode("producer")
