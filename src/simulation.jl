# ------------------------------------------------ #
# Simulation results structure
# ------------------------------------------------ #
"""Results of a simulation run.

`SimulationResults` stores the time series produced by [`run_simulation`](@ref).

```julia
struct SimulationResults
        time::Union{Vector{Float64}, Vector{DateTime}}
        mass_flow_load::Matrix{Float64}
        mass_flow_producer::Vector{Float64}
        T_load_in::Matrix{Float64}
        T_load_out::Matrix{Float64}
        T_producer_in::Vector{Float64}
        T_producer_out::Vector{Float64}
        power_load::Matrix{Float64}
        power_producer::Vector{Float64}
        load_labels::Dict{String, Int}
end
```

# Fields
- `time`: simulation time vector.
    - `Vector{Float64}`: time in seconds.
    - `Vector{DateTime}`: absolute timestamps.
- `mass_flow_load`: load mass flows in kg/s. Size `N × nloads`.
- `mass_flow_producer`: producer mass flow in kg/s. Length `N`.
- `T_load_in`: temperature entering each load (supply side) in °C. Size `N × nloads`.
- `T_load_out`: temperature leaving each load (return side) in °C. Size `N × nloads`.
- `T_producer_in`: return temperature entering the producer in °C. Length `N`.
- `T_producer_out`: supply temperature leaving the producer in °C. Length `N`.
- `power_load`: load power consumption in kW. Size `N × nloads`.
- `power_producer`: producer power output in MW (computed from mass flow and ΔT). Length `N-1`.
- `load_labels`: mapping from load label to column index used in the `*_load` matrices.

# Indexing
Convenience accessors are provided:

- `sr[:time]` returns the time vector.
- `sr[:load_labels]` returns the load labels.
- `sr[:load_labels_dict]` returns the label→column dictionary.
- `sr["L1", :T_load_in]` returns the time series for that load L1 (a vector).

# Notes
- All matrices are organized as `(time step, load index)`.
- `power_producer` has length `N-1` because the producer heats the water that entered in in *previous!* time step.
"""
struct SimulationResults
    time::Union{Vector{Float64}, Vector{DateTime}}  # time vector
    mass_flow_load::Matrix{Float64}         # mass flows at load nodes (rows: time steps, columns: load nodes)
    mass_flow_producer::Vector{Float64}     # mass flow at producer node
    T_load_in::Matrix{Float64}              # temperatures at load nodes entering (rows: time steps, columns: load nodes)
    T_load_out::Matrix{Float64}             # temperatures at load nodes exiting (rows: time steps, columns: load nodes)
    T_producer_in::Vector{Float64}          # input temperature entering producer node (after backward simulation step)
    T_producer_out::Vector{Float64}         # output temperature exiting producer node (before forward simulation step)
    power_load::Matrix{Float64}             # (kW) power consumption at load nodes (rows: time steps, columns: load nodes)
    power_producer::Vector{Float64}         # (MW) power output at producer node
    load_labels::Dict{String, Int}          # labels of load nodes corresponding to columns
end
# constructor with keyword arguments for better readability
function SimulationResults(;time, mass_flow_load, mass_flow_producer, T_load_in, T_load_out, T_producer_in, T_producer_out, power_load, power_producer, load_labels_dict)
    return SimulationResults(time, mass_flow_load, mass_flow_producer, T_load_in, T_load_out, T_producer_in, T_producer_out, power_load, power_producer, load_labels_dict)
end

# ------------------------------------------------ #
# Producer (power plant) output
# ------------------------------------------------ #
"""Control input returned by a simulation policy.

`ProducerOutput` represents the producer setpoints for one time step.

```julia
struct ProducerOutput
    mass_flow::Float64
    temperature::Float64
end
```

# Fields
- `mass_flow`: total mass flow injected into the network in kg/s.
- `temperature`: producer outlet (supply) temperature in °C.

# Usage
The `policy` passed to [`run_simulation`](@ref) must return a `ProducerOutput`:

```julia
function policy(t, Tₐ, T_back)
    return ProducerOutput(mass_flow=15.0, temperature=90.0)
end
```
"""
struct ProducerOutput
    mass_flow::Float64
    temperature::Float64
end
ProducerOutput(;mass_flow, temperature) = ProducerOutput(mass_flow, temperature)

function Base.getindex(sr::SimulationResults, label::String, s::Symbol)
    # usage: sr["M2_VS_1", :mass_flow] to get mass flow time series for load node "M2_VS_1"
    options = setdiff(fieldnames(SimulationResults), [:load_labels])
    if(s ∉ options && s != :load_labels_dict )
        error("$(s): Invalid symbol for SimulationResults getindex. Valid symbols are: $(options)")
    end
    idx = sr.load_labels[label]  # get column index for the load label
    if isnothing(idx)
        error("Load label $label not found in SimulationResults.")
    end
    if occursin("load", String(s))
        return getproperty(sr, s)[:, idx]
    else
        return getproperty(sr, s)
    end
end
function Base.getindex(sr::SimulationResults, s::Symbol)
    # usage example: sr[:time] to get time vector
    options = (fieldnames(SimulationResults)..., :load_labels_dict)
    if(s ∉ options)
        error("$(s): Invalid symbol for SimulationResults getindex. Valid symbols are: $(options)")
    end
    if(s==:load_labels)
        return keys(sr.load_labels)
    elseif (s==:load_labels_dict)
        return sr.load_labels
    end
    return getproperty(sr, s)
end


Base.length(sr::SimulationResults) = length(sr.time)
Base.show(io::IO, sr::SimulationResults) = print(io, "SimulationResults with $(length(sr)) time steps and $(length(sr.load_labels)) load nodes.")

# ------------------------------------------------ #
# SIMULATION FUNCTIONS
# ------------------------------------------------ #

"""Run a quasi-dynamic simulation of a district heating network.

```julia
run_simulation(network, sim_time, policy; T0_f=60.0, T0_b=25.0, ambient_temperature=nothing)
```

This is the main entry point for time stepping.

REPEAT for N time steps:
1. computes a steady-state hydraulic solution (mass flow distribution),
2. advances thermal dynamics using the plug-flow method:
     - forward/supply advection producer → loads,
     - heat consumption at loads,
     - backward/return advection loads → producer,
     - heat losses to ambient.

See [Plug method](@ref) for the underlying model.


# Output
- `SimulationResults` struct containing time series of temperatures, flows, and powers for all nodes and edges.

# Arguments
- `network::Network`: prepared network (producer/load nodes identified, pipes attached).
- `sim_time`: equally spaced time vector.
    - `Vector{Float64}`: time in seconds.
    - `Vector{DateTime}`: timestamps (Δt is interpreted in seconds).
- `policy::Function`: callback returning `ProducerOutput`.
    - Signature: `policy(t, Tₐ, T_back)::ProducerOutput`
    - `Tₐ` is ambient temperature at time `t` or `nothing`.
    - `T_back` is the return temperature entering the producer (in previous time step `k-1`).

# Keyword Arguments
- `T0_f`: initial temperature forward part of the network (producer → loads) (°C).
- `T0_b`: initial temperature in backward part of the network (loads → producer) (°C).
- `ambient_temperature`: optional `Vector{Float64}` of ambient (outdoor/atmospheric) temperatures (°C), length must match `sim_time`.

# Returns
- `SimulationResults`: time series of temperatures, flows, and powers.

# Notes
- The network structure is validated once at the start via `check_network!`.
- Time steps must be equally spaced.
"""
function run_simulation(network::Network, sim_time::Union{Vector{Float64}, Vector{DateTime}}, policy::Function; 
                        T0_f::Float64=60.0, T0_b::Float64=25.0, ambient_temperature::Union{Vector{Float64}, Nothing}=nothing)::SimulationResults
    # simulate dynamics of the DH Network
    # results are y(t) = f(x(t), u(t)), where u(t) are inputs: (mass_flow, input_temperature) for source node
    # x0 is initial state of the network ... in default variant we fill the pipes with water of 25 °C
    # inputs are (mass_flow_input, temp_input) vectors of length N (number of time steps)
    # policy(t, Tₐ, T_back)::ProducerOutput
    #   t is current time (in seconds or DateTime)
    #   Tₐ is ambient temperature at time t or nothing if ambient temperature is not provided
    #   T_back is the temperature of water returning to producer at time t (if t > 1, otherwise use initial T0_b)

    # check network structure (one producer, connected, acyclic,...), also updates neighbor dicts if they are not built yet
    check_network!(network::Network)

    # check time steps are equally spaced
    if(eltype(sim_time) <: Float64)
        dt = diff(sim_time)
    elseif (eltype(sim_time) <: DateTime)
        dt = diff(sim_time) .|> Dates.Second .|> Dates.value # compute time steps in seconds and convert to Int
    else
        error("Unsupported time vector element type: $(eltype(sim_time)). Expected Float64 or DateTime.")
    end
    if !all(isapprox.(dt, dt[1], atol=1e-10))
        error("Time steps are not equally spaced!")
    end
    Δt = float(dt[1])

    # ckeck policy function output for the first time step
    try
        Tₐ = isnothing(ambient_temperature) ? nothing : ambient_temperature[1]
        if policy(sim_time[1], Tₐ, T0_b) isa ProducerOutput
            # ok
        else
            error("Policy function must return a ProducerOutput struct!")
        end
    catch e
        error("Error when calling policy function for the first time step: ", e)
    end

    # fill pipes with initial temperatures
    fill_pipes_with_initial_temperature!(network, T0_f, T0_b)

    N = length(sim_time)
    num_loads = length(network.load_labels)

    # store results: each column corresponds to a load node, each row to a time step
    results_mass_flow_producer = Vector{Float64}(undef, N)
    results_mass_flow_load = Matrix{Float64}(undef, N, num_loads)
    results_temperature_producer_out = Vector{Float64}(undef, N)
    results_temperature_producer_in = Vector{Float64}(undef, N)
    results_temperature_load_in = Matrix{Float64}(undef, N, num_loads)
    results_temperature_load_out = Matrix{Float64}(undef, N, num_loads)
    results_power_consumption = Matrix{Float64}(undef, N, num_loads)
    
    # mapping from load node label to column index in results matrices...
    load_labels_cols = Dict(label => i for (i, label) in enumerate(network.load_labels)) 

    
    for i in 1:N
        # GET INPUTS FOR THIS TIME STEP
        Tₐ = isnothing(ambient_temperature) ? nothing : ambient_temperature[i]
        T_back = i > 1 ? results_temperature_producer_in[i-1] : T0_b 
        input = policy(sim_time[i], Tₐ, T_back)

        results_temperature_producer_out[i] = input.temperature
        results_mass_flow_producer[i] = input.mass_flow


        # FORWARD SIMULATION STEP
        steady_state_hydronynamics!(network, input.mass_flow)
        Tₐ = isnothing(ambient_temperature) ? nothing : ambient_temperature[i]
        output_plugs = time_step_thermal_dynamics_forward!(network, Δt, input.temperature, Tₐ)
        
        # log output values
        for (load_label, plug) in output_plugs
            col_idx = load_labels_cols[load_label]
            results_temperature_load_in[i, col_idx] = plug.T
            results_mass_flow_load[i, col_idx] = plug.m / Δt
        end

        # POWER CONSUMPTION STEP
        Tₐ = isnothing(ambient_temperature) ? 15.0 : ambient_temperature[i] # use 15 °C as default ambient temperature
        for load_label in keys(output_plugs)
            col_idx = load_labels_cols[load_label]
            # power depends on outdoor temperature
            P = power_consumption(network[load_label], Tₐ)
            results_power_consumption[i, col_idx] = P / 1000.0 # log power consumption in kW
            
            # cool plug according to power consumed before feeding it back to return flow
            consume_power!(output_plugs[load_label], P, Δt)
            results_temperature_load_out[i, col_idx] = output_plugs[load_label].T
        end

        # BACKWARD SIMULATION STEP
        incoming_plug = time_step_thermal_dynamics_backward!(network, Δt, output_plugs, Tₐ)
        results_temperature_producer_in[i] = incoming_plug.T
    end

    # compute producer power output in MW based on mass flow and temperature difference between producer input and output
    # P[k] = ̇m[k] * c * (T_out[k] - T_in[k-1]) / 1_000_000.0 to convert from W to MW
    power_producer = @. (results_temperature_producer_out[2:end] - results_temperature_producer_in[1:end-1]) * results_mass_flow_producer[1:end-1] * WATER_SPECIFIC_HEAT / 1_000_000.0 # in MW

    return SimulationResults(time=sim_time,
                            mass_flow_load=results_mass_flow_load,
                            mass_flow_producer=results_mass_flow_producer,
                            T_load_in=results_temperature_load_in,
                            T_load_out=results_temperature_load_out,
                            T_producer_in=results_temperature_producer_in,
                            T_producer_out=results_temperature_producer_out, 
                            power_load=results_power_consumption,
                            power_producer=power_producer, 
                            load_labels_dict=load_labels_cols)
end

# ------------------------------------------------- #
# Simulation helper methods
# ------------------------------------------------- #

function set_relative_mass_flows!(nw::Network)
    # iterate over nodes from leaves to root and set on each edge relative mass flow coefficient m_rel
    # we will do iterative post-order DFS traversal

    root = nw.producer_label
    if root === nothing
        error("Producer node is not set in the network (producer_label is nothing).")
    end
    root = root::String
    stack = Vector{Tuple{String, Bool}}()  # stack of (node, visited)
    push!(stack, (root, false))
    while !isempty(stack)
        node, visited = pop!(stack)

        if visited
            # SECOND VISIT -> process node
            if nw[node] isa ProducerNode
                continue # skip root node
            end

            parent_node = inneighbors(nw, node)[1] # assume there is only one parent
            if(outdegree(nw, node) == 0) # leaf node
                @assert nw[node] isa LoadNode
                nw[parent_node, node].m_rel = nw[node].m_rel # copy m_rel from LoadNode to the edge leading to it
                    
            elseif nw[node] isa JunctionNode # visiting for second time, so all childer have been processed
                m_rel_sum = sum(nw[node, child].m_rel for child in outneighbors(nw, node))::Float64
                nw[parent_node, node].m_rel = m_rel_sum # set m_rel on the edge leading to this junction
            end

        else
            # FIRST VISIT -> defer node
            push!(stack, (node, true))

            # push children (unvisited)
            for child in outneighbors(nw, node)
                push!(stack, (child, false))
            end
        end
    end
end

function set_absolute_mass_flows!(nw::Network, mass_flow_source::Float64)
    # set absolute mass flows on edges based on relative mass flow coefficients m_rel and source mass flow
    # when entering junction, mass flow is split according to m_rel coefficients
    # do it by BFS from root to leaves

    root = nw.producer_label
    if root === nothing
        error("Producer node is not set in the network (producer_label is nothing).")
    end
    root = root::String
    nw[root].common.mass_flow = mass_flow_source
    queue = [root]

    while !isempty(queue)
        node = popfirst!(queue)
        if(outdegree(nw, node) == 0) # leaf node
            continue
        end

        sum_rel_flow = sum(nw[node, child].m_rel for child in outneighbors(nw, node))
        if sum_rel_flow == 0.0
            warn("Node $(nw[node].common.info) has zero total relative mass flow to its children. Skipping mass flow assignment for its outgoing edges.")
            continue # no flow through this node
        end
        for child in outneighbors(nw, node)
            edge = nw[node, child]
            edge.mass_flow = nw[node].common.mass_flow * edge.m_rel/sum_rel_flow
            nw[child].common.mass_flow = edge.mass_flow
            push!(queue, child)
        end
    end
end

function steady_state_hydronynamics!(nw::Network, mass_flow_source::Float64)
    # compute steady-state hydrodynamics of the network
    # assuming that mass flows are already set on edges and nodes

    set_relative_mass_flows!(nw)
    set_absolute_mass_flows!(nw, mass_flow_source)
end

function fill_pipes_with_initial_temperature!(nw::Network, temperature_f::Float64, temperature_b::Float64)
    # fill all pipes in the network with water at initial temperature, forward and backward have different T0
    for e in edges(nw.mg)
        edge = nw[e.src, e.dst]
        if edge isa InsulatedPipe
            L = pipe_length(edge)  # length in m
            d = inner_diameter(edge)  # inner diameter in m
            V = π * (d/2)^2 * L  # volume in m^3
            m_total = WATER_DENSITY * V  # total mass in kg
            empty!(edge.plugs_f)
            empty!(edge.plugs_b)
            push!(edge.plugs_f, Plug(temperature_f, m_total))
            push!(edge.plugs_b, Plug(temperature_b, m_total))
        end
    end
end

function time_step_thermal_dynamics_forward!(nw::Network, Δt::Float64, temperature_source::Float64, ambient_temperature::Union{Float64, Nothing}=nothing)::Dict{String, Plug}
    # simulate thermal dynamics of the network for one time step Δt
    # update temperatures in plugs in all pipes based on heat losses and advection

    # create plugs entering the network at the source node
    root = nw.producer_label
    if root === nothing
        error("Producer node is not set in the network (producer_label is nothing).")
    end
    root = root::String
    source_edges = [nw[root, n] for n in outneighbors(nw, root)]
    plug_masses = [edge.mass_flow * Δt for edge in source_edges] # integrate mass flow from source edge to get mass in kg
    for (source_edge, plug_mass) in zip(source_edges, plug_masses)
        new_plug = Plug(temperature_source, plug_mass) # plug entering the first edge
        push_in_water_plugs_forward!(source_edge, [new_plug])
    end

    output_plugs = Dict{String, Plug}()  # dictionary of plugs exiting each leaf node

    # in BFS manner, process all edges
    queue = [root]
    for child in outneighbors(nw, root)
        push!(queue, child)
    end

    while !isempty(queue)
        node = popfirst!(queue)
        if indegree(nw, node) == 0 # root node
            # already processed above
            continue
        end
        parent_edge = nw[inneighbors(nw, node)[1], node] # assume there is only one parent edge
        
        # collect plugs exiting from parent edge
        plugs = collect_exiting_water_plugs!(parent_edge.plugs_f, parent_edge.mass_flow, Δt)
        
        if(outdegree(nw, node) == 0) # leaf node
            # log exiting plugs for this load node (average them weighting by mass)
            output_plugs[node] = combine_plugs(plugs)
            continue
        else
            # push child nodes to queue
            for child in outneighbors(nw, node)
                push!(queue, child)
            end
        end

        # devide each plug into child edges according to mass flow in each edge
        children = outneighbors(nw, node)
        total_mass_flow = sum(nw[node, child].mass_flow for child in children)
        
        if(total_mass_flow == 0.0 && !isempty(plugs))
            error("Node $(nw[node].common.info) has zero total mass flow to its children. Cannot distribute plugs!")
        end
        for child in children
            edge = nw[node, child]
            mass_flow_edge = edge.mass_flow
            next_pipe_plugs = Vector{Plug}()
            for p in plugs
                next_mass = p.m * mass_flow_edge / total_mass_flow
                if next_mass > 0.0
                    next_plug = Plug(p.T, next_mass)
                    push!(next_pipe_plugs, next_plug)
                end
            end
            # push the smaller plugs into the child edge
            push_in_water_plugs_forward!(edge, next_pipe_plugs)
        end
    end

    # cool remaining plugs in the network due to heat loss
    if ambient_temperature !== nothing
        heat_loss_forward!(nw, ambient_temperature, Δt)
    end

    return output_plugs
end

function time_step_thermal_dynamics_backward!(nw::Network, Δt::Float64, incoming_plugs::Dict{String, Plug}, ambient_temperature::Union{Float64, Nothing}=nothing)::Union{Plug, Nothing}
    # simulate thermal dynamics of the network for one time step Δt
    # fill back in the network the cooled plugs from load nodes and propagate up to the producer node
    # update temperatures in plugs in all pipes based on heat losses and advection
    output_plug = nothing

    # iterate over nodes from leaves to root using iterative post-order DFS traversal

    root = nw.producer_label
    if root === nothing
        error("Producer node is not set in the network (producer_label is nothing).")
    end
    root = root::String
    stack = Vector{Tuple{String, Bool}}()  # stack of (node, visited)
    push!(stack, (root, false))
    while !isempty(stack)
        node, visited = pop!(stack)

        if visited
            # SECOND VISIT -> process node
            if (outdegree(nw,node)) == 0 # leaf node
                @assert nw[node] isa LoadNode
                parent = inneighbors(nw, node)[1] # assume there is only one parent
                parent_edge = nw[parent, node]
                push_in_water_plugs_backward!(parent_edge, [incoming_plugs[node]]) # push the cooled plug back into the parent edge
                    
            else # visiting for second time, so all childer have been processed
                children_plug_vectors = Vector{Vector{Plug}}() # collect plugs from child edges to merge them
                # collect plugs exiting from childs edges
                for child in outneighbors(nw, node)
                    edge = nw[node, child]
                    plugs = collect_exiting_water_plugs!(edge.plugs_b, edge.mass_flow, Δt)
                    push!(children_plug_vectors, plugs)
                end
                # combine plugs from all child edges
                merged = merge_water_plug_vectors!(children_plug_vectors)
                # push merged plugs back to parent edge
                if node != root
                    parent = inneighbors(nw, node)[1] # assume there is only one parent
                    parent_edge = nw[parent, node]
                    push_in_water_plugs_backward!(parent_edge, merged)
                else
                    # at the root, we can return the combined plug as output of backward simulation
                    output_plug = combine_plugs(merged)
                    break
                end
            end

        else
            # FIRST VISIT -> defer node
            push!(stack, (node, true))

            # push children (unvisited)
            for child in outneighbors(nw, node)
                push!(stack, (child, false))
            end
        end
    end

    # cool remaining plugs in the network due to heat loss
    if !isnothing(ambient_temperature)
        heat_loss_backward!(nw, ambient_temperature, Δt)
    end

    return output_plug
end

function collect_exiting_water_plugs!(plugs::Vector{Plug}, mass_flow::Float64, Δt::Float64)::Vector{Plug}
    # collect plugs that exit the pipe during time step Δt
    exited_plugs = Vector{Plug}()
    mass_exited = mass_flow * Δt # theoretical amount that should exit
    if(mass_exited <= 0.0)
        return exited_plugs  # no mass exiting
    end
    mass_accumulated = 0.0
    while !isempty(plugs) && mass_accumulated < mass_exited
        p = popfirst!(plugs)
        if mass_accumulated + p.m > mass_exited
            # only part of the plug exits
            remaining_mass = mass_exited - mass_accumulated
            exiting_plug = Plug(p.T, remaining_mass)
            push!(exited_plugs, exiting_plug)
            p.m -= remaining_mass
            pushfirst!(plugs, p) # put the remaining part of the plug back to the front of the queue
            mass_accumulated = mass_exited
            break
        else # entire plug exits
            push!(exited_plugs, p)
            mass_accumulated += p.m
        end
    end
    if isempty(plugs) && !isapprox(mass_accumulated, mass_exited; atol=1e-6) # this should not happen
        error("Not enough mass in the pipe to exit! Mass exited: $mass_accumulated, expected: $mass_exited. This indicates a problem with the simulation, possibly due to numerical errors or incorrect mass flow assignment.")
    end
    return exited_plugs
end
function push_in_water_plugs_forward!(e::InsulatedPipe, plugs::Vector{Plug})
    # push incoming plugs into the pipe's plug queue
    for p in plugs
        push!(e.plugs_f, p)
    end
    merge_same_temperature_plugs!(e.plugs_f)
end
function push_in_water_plugs_backward!(e::InsulatedPipe, plugs::Vector{Plug})
    # push incoming plugs into the pipe's plug queue
    for p in plugs
        push!(e.plugs_b, p)
    end
    merge_same_temperature_plugs!(e.plugs_b)
end

function merge_same_temperature_plugs!(plugs::Vector{Plug}; tol::Float64=1e-3)
    # merge consecutive plugs with same temperature (within tolerance tol) to reduce number of plugs
    if isempty(plugs)
        return
    end
    merged_plugs = Vector{Plug}()
    current_plug = popfirst!(plugs)
    while !isempty(plugs)
        next_plug = popfirst!(plugs)
        if abs(current_plug.T - next_plug.T) < tol
            # merge plugs
            current_plug = combine_plugs([current_plug, next_plug])
        else
            push!(merged_plugs, current_plug)
            current_plug = next_plug
        end
    end
    push!(merged_plugs, current_plug)
    empty!(plugs)
    append!(plugs, merged_plugs) # replace original plugs with merged ones
end

function combine_plugs(plugs::Vector{Plug})::Plug
    # combine multiple plugs into one by mass-weighted average of temperature
    total_mass = sum(p.m for p in plugs)
    if total_mass == 0.0
        return Plug(25.0, 0.0) # default temperature for zero mass
    end
    avg_temp = sum(p.T * p.m for p in plugs) / total_mass
    return Plug(avg_temp, total_mass)
end

function set_load_params!(nw::Network, load_label::String, m_rel::Float64)
    # set parameters of a load node, for now just relative mass flow coefficient
    @assert nw[load_label] isa LoadNode
    nw[load_label].m_rel = m_rel
end

function set_load_params!(nw::Network, load_params::Dict{String, Float64})
    # set parameters of multiple load nodes
    for (load_label, m_rel) in load_params
        set_load_params!(nw, load_label, m_rel)
    end
end

function heat_loss_forward!(e::InsulatedPipe, T_a::Float64, Δt::Float64)
    # compute heat loss in the pipe and cool the plugs accordingly
    ρ = WATER_DENSITY  # density in kg/m^3
    cₚ = WATER_SPECIFIC_HEAT  # specific heat capacity in J/(kg·K)
    for p in e.plugs_f
        A = π * (inner_diameter(e)/2)^2  # surface area in m^2
        T_next = T_a + (p.T - T_a) * exp(- (Δt) / (ρ * cₚ * A * heat_resistance_forward(e)))
        p.T = T_next
    end
end

function heat_loss_forward!(nw::Network, ambient_temperature::Float64, Δt::Float64)
    # compute heat loss in all pipes of the network and cool the plugs accordingly
    for e in edges(nw.mg)
        edge = nw[e.src, e.dst]
        if edge isa InsulatedPipe
            heat_loss_forward!(edge, ambient_temperature, Δt)
        end
    end
end

function heat_loss_backward!(e::InsulatedPipe, T_a::Float64, Δt::Float64)
    # compute heat loss in the pipe and cool the plugs accordingly
    ρ = WATER_DENSITY  # density in kg/m^3
    cₚ = WATER_SPECIFIC_HEAT  # specific heat capacity in J/(kg·K)
    for p in e.plugs_b
        A = π * (inner_diameter(e)/2)^2  # surface area in m^2
        T_next = T_a + (p.T - T_a) * exp(- (Δt) / (ρ * cₚ * A * heat_resistance_backward(e)))
        p.T = T_next
    end
end

function heat_loss_backward!(nw::Network, ambient_temperature::Float64, Δt::Float64)
    # compute heat loss in all pipes of the network and cool the plugs accordingly
    for e in edges(nw.mg)
        edge = nw[e.src, e.dst]
        if edge isa InsulatedPipe
            heat_loss_backward!(edge, ambient_temperature, Δt)
        end
    end
end

function power_consumption(node::LoadNode, Tₐ::Float64)::Float64
    # compute power consumption of a load node based on outdoor temperature and load coefficients
    # power consumtion is polynomial function of outdoor temperature: P = a + b*Tₐ + c*Tₐ^2 + ...
    T₀ = -node.load[2] / (2*node.load[3]) # minimal power consumption at this outdoor temperature
    if(Tₐ > T₀) # function must be decresing
        return node.load[1] + node.load[2]*T₀ + node.load[3]*T₀^2
    end

    P = 0.0
    for i in 1:length(node.load)
        P += node.load[i] * (Tₐ^(i-1))
    end
    return P * 1000.0 # convert from kW to W
end

function consume_power!(p::Plug, power::Float64, Δt::Float64)
    # compute new temperature of the plug after consuming power for time step Δt
    # energy consumed is E = P * Δt, which reduces the thermal energy of the plug: m*cₚ*ΔT = E
    ρ = WATER_DENSITY  # density in kg/m^3
    cₚ = WATER_SPECIFIC_HEAT  # specific heat capacity in J/(kg·K)
    ΔT = power * Δt / (p.m * cₚ)  # temperature drop in K
    if(p.T - ΔT < MINIMAL_RETURN_TEMPERATURE)
        @warn "Return temperature at load too low, limiting to minimal return temperature of $MINIMAL_RETURN_TEMPERATURE °C to avoid unphysical results."
        ΔT = p.T - MINIMAL_RETURN_TEMPERATURE # limit temperature drop to avoid unphysical results
    end
    p.T -= ΔT
end

function merge_water_plug_vectors!(plug_vectors::Vector{Vector{Plug}})::Vector{Plug}
    # merge multiple vectors of plugs into one vector by combining plugs with same temperature
    # first we need to find the times when plugs change ... parametrize by t ∈ [0,1]

    if length(plug_vectors) == 1
        merge_same_temperature_plugs!(plug_vectors[1]; tol=1e-2)
        return plug_vectors[1] # no need to merge if there is only one vector
    end


    change_points = Set{Float64}()
    total_masses = [sum(p.m for p in plugs) for plugs in plug_vectors]
    for (i, plugs) in enumerate(plug_vectors)
        mass_accumulated = 0.0
        for p in plugs
            mass_accumulated += p.m
            push!(change_points, mass_accumulated / total_masses[i])
        end
    end

    change_points = sort(collect(change_points))


    # now combine plugs between change points by mass-weighted average of temperature
    merged_plugs = Vector{Plug}()
    t_prev = 0.0
    for t in sort(collect(change_points))
        plugs_at_t = Vector{Plug}()
        for i in eachindex(plug_vectors)
            target_mass = (t-t_prev) * total_masses[i]
            p = popfirst!(plug_vectors[i])
            if !isapprox(p.m, target_mass; atol=1e-6)
                # only part of the plug is in this interval
                # if remaining part is very small, throw away
                pushfirst!(plug_vectors[i], Plug(p.T, p.m - target_mass)) # put the remaining part back to the front of the queue
                p = Plug(p.T, target_mass) # take only the part of the plug that is in this interval
            end
            push!(plugs_at_t, p)
        end
        t_prev = t
        push!(merged_plugs, combine_plugs(plugs_at_t)) # combine plugs at this change point into one plug
    end

    # there should no plugs remain in the original vectors
    @assert all(isempty, plug_vectors)

    merge_same_temperature_plugs!(merged_plugs; tol=1e-2) # simplify the result by merging plugs with almost the same temperature
    
    @assert isapprox(sum(p.m for p in merged_plugs), sum(total_masses); atol=1e-6) # check that total mass is conserved
    return merged_plugs
end