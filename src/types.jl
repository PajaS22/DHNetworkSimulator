# ------------------------------------------------ #
# TYPE DEFINITIONS
# ------------------------------------------------ #

# abstract types for DH Network Nodes and Edges
abstract type dhNodeType end
abstract type dhEdgeType end

# Common Node Data (shared by all node types)
"Common data for all node types in the DH network (e.g., junctions, loads, producers)"
mutable struct dhNodeCommon
    info::String
    position::Union{Missing, Tuple{Float64, Float64}}  # (x, y) coordinates
    mass_flow::Union{Missing, Float64}      # [kg/s]
end

# Junction Node: connection point with no heat consumption/production
"DH network node representing a junction, where pipes can connect but no heat is produced or consumed."
struct dhJunctionNode <: dhNodeType
    common::dhNodeCommon
end

# Load Node: heat consumer
"DH network node representing a load, which consumes heat. The load field defines quadratic function of power consumption P(Tₐ) = p₀ + p₁*Tₐ + p₂*Tₐ²."
mutable struct dhLoadNode <: dhNodeType
    common::dhNodeCommon
    load::Union{Missing, NTuple{3, Float64}}     # Heat load function in kW, P(Tₐ) = p₀ + p₁*Tₐ + p₂*Tₐ²
    m_rel::Union{Missing, Float64}               # Relative mass flow coefficient (for branching nodes)
end

# Producer Node: heat producer
"DH network node representing a producer, which heats water. There may be only one producer in the network, and it is identified by its label (producer_label field in dhNetwork)."
struct dhProducerNode <: dhNodeType
    common::dhNodeCommon
end

# empty Node: just for initialization purposes
"An empty node type used for initialization and placeholder purposes in the DH network."
struct EmptyNode <: dhNodeType end

# ------------------------------------------------ #
# PIPE MODEL
# ------------------------------------------------ #
# we will model pipe using "plug" method

mutable struct Plug
    T::Float64  # Temperature at the plug [°C]
    m::Float64  # mass of the plug [kg]
end

struct PipeParams                       # unchanging physical parameters of the pipe
    length::Float64                     # Length of the pipe [m]
    inner_diameter::Float64             # Inner diameter [m]
    heat_resistance_forward::Float64    # Thermal resistance [m*K/W]
    heat_resistance_backward::Float64   # Thermal resistance [m*K/W]
end

# Pipe Edge: represents a pipe in the network
mutable struct dhPipeEdge <: dhEdgeType
    info::String
    physical_params::PipeParams
    mass_flow::Union{Missing, Float64}  # Mass flow in [kg/s]
    m_rel::Union{Missing, Float64}      # Relative mass flow coefficient (for branching pipes)
    plugs_f::Vector{Plug}               # Queue of plugs in the pipe (forward direction)
    plugs_b::Vector{Plug}               # Queue of plugs in the pipe (backward direction)
end

struct EmptyEdge <: dhEdgeType end

# ------------------------------------------------- #
# DH NETWORK TYPE
# ------------------------------------------------- #

mutable struct dhNetwork{T<:Integer} <: AbstractGraph{T}
    mg::MetaGraph                               # MetaGraph from MetaGraphs.jl, it contains the network structure and node and edge data
    producer_label::Union{Nothing, String}      # Label of the producer node
    load_labels::Set{String}                    # Labels of the consumer nodes
end

# ------------------------------------------------- #
# CONSTRUCTORS
# ------------------------------------------------- #

# default constructor for empty network, more about construction in network_creation.jl
function dhNetwork()
    mg = MetaGraph(
        DiGraph();  # underlying graph structure
        label_type=String,
        vertex_data_type=dhNodeType,
        edge_data_type=dhEdgeType
    )
    return dhNetwork{Int}(mg, nothing, Set{String}())
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
function InsulatedPipe(info::String, params::PipeParams)::dhPipeEdge
    mass_flow = missing         # [kg/s]
    m_rel = missing             # Relative mass flow coefficient
    plugs_f = Vector{Plug}()    # Initialize empty queue of plugs (forward direction)
    plugs_b = Vector{Plug}()    # Initialize empty queue of plugs (backward direction)
    return dhPipeEdge(info, params, mass_flow, m_rel, plugs_f, plugs_b)
end
InsulatedPipe(params::PipeParams) = InsulatedPipe("pipe", params)
InsulatedPipe(length::Real) = InsulatedPipe("pipe"; length=float(length))

# ------------------------------------------------- #
# NODE CONSTRUCTORS
# ------------------------------------------------- #

# dhNodeCommon constructors
dhNodeCommon(info::String) = dhNodeCommon(info, missing, missing)
dhNodeCommon(info::String, position::Tuple{Float64, Float64}) = dhNodeCommon(info, position, missing)

# JUNCTION NODE CONSTRUCTORS
JunctionNode(info::String) = dhJunctionNode(dhNodeCommon(info))
JunctionNode(info::String, position::Tuple{Float64, Float64}) = dhJunctionNode(dhNodeCommon(info, position))
JunctionNode(position::Tuple{Float64, Float64}) = JunctionNode("junction", position)
JunctionNode() = JunctionNode("junction")

# LOAD NODE CONSTRUCTORS
const DEFAULT_LOAD = (540.0, -36.0, 0.6)
LoadNode(info::String; load=DEFAULT_LOAD) = dhLoadNode(dhNodeCommon(info), load, missing)
LoadNode(info::String, position::Tuple{Float64, Float64}; load=DEFAULT_LOAD) = dhLoadNode(dhNodeCommon(info, position), load, missing)
LoadNode(info::String, position::Tuple{Float64, Float64}, load::NTuple{3, Float64}) = dhLoadNode(dhNodeCommon(info, position), load, missing)
LoadNode(info::String, position::Tuple{Float64, Float64}, m_rel::Float64; load=DEFAULT_LOAD) = dhLoadNode(dhNodeCommon(info, position), load, m_rel)
LoadNode(position::Tuple{Float64, Float64}; load=DEFAULT_LOAD) = LoadNode("load", position; load=load)
LoadNode(; load=DEFAULT_LOAD) = LoadNode("load"; load=load)

# PRODUCER NODE CONSTRUCTORS
ProducerNode(info::String) = dhProducerNode(dhNodeCommon(info))
ProducerNode(info::String, position::Tuple{Float64, Float64}) = dhProducerNode(dhNodeCommon(info, position))
ProducerNode(position::Tuple{Float64, Float64}) = ProducerNode("producer", position)
ProducerNode() = ProducerNode("producer")
