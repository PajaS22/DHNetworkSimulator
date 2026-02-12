# ------------------------------------------------ #
# TYPE DEFINITIONS
# ------------------------------------------------ #


"Abstract type for all node types in the DH network (e.g., junctions, loads, producers)."
abstract type NodeType end

"Abstract type for all edge types in the DH network (only one pipe for now)."
abstract type EdgeType end


"""Common data for all node types in the DH network (e.g., junctions, loads, producers)

    mutable struct NodeCommon
        info::String
        position::Union{Missing, Tuple{Float64, Float64}}   # (x, y) coordinates
        mass_flow::Union{Missing, Float64}                  # [kg/s]
    end

    # Constructors
    - `NodeCommon(info::String)`
    - `NodeCommon(info::String, position::Tuple{Float64, Float64})`
"""
mutable struct NodeCommon
    info::String
    position::Union{Missing, Tuple{Float64, Float64}}   # (x, y) coordinates
    mass_flow::Union{Missing, Float64}                  # [kg/s]
end


"""DH network node representing a junction, where pipes can connect but no heat is produced or consumed.

    struct JunctionNode <: NodeType
        common::NodeCommon
    end

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
"""DH network node representing a load, which consumes heat.

The load field defines quadratic function of power consumption P(Tₐ) = p₀ + p₁*Tₐ + p₂*Tₐ².

    mutable struct LoadNode <: NodeType
        common::NodeCommon
        load::Union{Missing, NTuple{3, Float64}}     # Heat load function in kW, P(Tₐ) = p₀ + p₁*Tₐ + p₂*Tₐ²
        m_rel::Union{Missing, Float64}               # Relative mass flow coefficient (for branching nodes)
    end

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
"""DH network node representing a producer, which heats water.

There may be only one producer in the network, and it is identified by its label (producer_label field in Network).

    struct ProducerNode <: NodeType
        common::NodeCommon
    end

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
"An empty node type used for initialization and placeholder purposes in the DH network."
struct EmptyNode <: NodeType end

# ------------------------------------------------ #
# PIPE MODEL
# ------------------------------------------------ #
# we will model pipe using "plug" method

"""
Single plug of water in the pipe, characterized by its temperature and mass.

    mutable struct Plug
        T::Float64  # Temperature at the plug [°C]
        m::Float64  # mass of the plug [kg]
    end
"""
mutable struct Plug
    T::Float64  # Temperature at the plug [°C]
    m::Float64  # mass of the plug [kg]
end

"""
Physical parameters of a pipe, which are constant during the simulation.

    struct PipeParams                       # unchanging physical parameters of the pipe
        length::Float64                     # Length of the pipe [m]
        inner_diameter::Float64             # Inner diameter [m]
        heat_resistance_forward::Float64    # Thermal resistance [m*K/W]
        heat_resistance_backward::Float64   # Thermal resistance [m*K/W]
    end

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
"""
DH network edge representing an insulated pipe, which transports water between nodes.

    mutable struct InsulatedPipe <: EdgeType
        info::String
        physical_params::PipeParams
        mass_flow::Union{Missing, Float64}  # Mass flow in [kg/s]
        m_rel::Union{Missing, Float64}      # Relative mass flow coefficient (for branching pipes)
        plugs_f::Vector{Plug}               # Queue of plugs in the pipe (forward direction)
        plugs_b::Vector{Plug}               # Queue of plugs in the pipe (backward direction)
    end

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

"Empty edge type used for initialization and placeholder purposes in the DH network."
struct EmptyEdge <: EdgeType end

# ------------------------------------------------- #
# DH NETWORK TYPE
# ------------------------------------------------- #

"Mappings to inneghbors and outneighbors for efficient access during simulation, stored in the Network struct."
mutable struct NeighborDicts
    outneighbors::Dict{String, Vector{String}}       # mapping nodes to outneighbors for efficient access in simulation
    inneighbors::Dict{String, Vector{String}}        # mapping nodes to inneighbors for efficient access in simulation
    need_rebuild::Bool                               # flag to indicate if neighbor dicts need to be rebuilt before simulation
end
NeighborDicts() = NeighborDicts(Dict{String, Vector{String}}(), Dict{String, Vector{String}}(), true)


"""
The main data structure representing a district heating network, which consists of a MetaGraph containing
nodes and edges with associated data, as well as labels for the producer and consumer nodes.

    mutable struct Network{T<:Integer} <: AbstractGraph{T}
        mg::MetaGraph                               # MetaGraph from MetaGraphs.jl, it contains the network 
                                                      structure and node and edge data
        producer_label::Union{Nothing, String}      # Label of the producer node
        load_labels::Set{String}                    # Labels of the consumer nodes
        neighbor_dicts::NeighborDicts               # Mappings to inneghbors and outneighbors for efficient access during simulation
    end

    # Constructors
    - `Network()`: Creates an empty DH network with no nodes or edges.
    - `Network(g::DiGraph)`: Creates a DH network from an existing directed graph structure.
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
