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
        T_load_out::Union{Matrix{Float64}, Nothing}
        T_producer_in::Union{Vector{Float64}, Nothing}
        T_producer_out::Vector{Float64}
        power_load::Union{Matrix{Float64}, Nothing}
        power_producer::Union{Vector{Float64}, Nothing}
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
- `T_load_out`: temperature leaving each load (return side) in °C. Size `N × nloads`. `Nothing` in forward-only mode.
- `T_producer_in`: return temperature entering the producer in °C. Length `N`. `Nothing` in forward-only mode.
- `T_producer_out`: supply temperature leaving the producer in °C. Length `N`.
- `power_load`: load power consumption in kW. Size `N × nloads`. `Nothing` in forward-only mode.
- `power_producer`: producer power output in MW (computed from mass flow and ΔT). Length `N-1`. `Nothing` in forward-only mode.
- `load_labels`: mapping from load label to column index used in the `*_load` matrices.

# Indexing
Convenience accessors are provided:

- `sr[:time]` returns the time vector.
- `sr[:load_labels]` returns the load labels.
- `sr[:load_labels_dict]` returns the label→column dictionary.
- `sr["L1", :T_load_in]` returns the time series for that load L1 (a vector).

# Notes
- All matrices are organized as `(time step, load index)`.
- `power_producer` has length `N-1` because the producer heats water that arrived from the *previous* time step — so there is one fewer value than time steps.
"""
struct SimulationResults
    time::Union{Vector{Float64}, Vector{DateTime}}  # time vector
    mass_flow_load::Matrix{Float64}         # mass flows at load nodes (rows: time steps, columns: load nodes)
    mass_flow_producer::Vector{Float64}     # mass flow at producer node
    T_load_in::Matrix{Float64}              # temperatures at load nodes entering (rows: time steps, columns: load nodes)
    T_load_out::Union{Matrix{Float64}, Nothing}     # temperatures at load nodes exiting (rows: time steps, columns: load nodes); Nothing in forward-only mode
    T_producer_in::Union{Vector{Float64}, Nothing}  # input temperature entering producer node (after backward simulation step); Nothing in forward-only mode
    T_producer_out::Vector{Float64}         # output temperature exiting producer node (before forward simulation step)
    power_load::Union{Matrix{Float64}, Nothing}     # (kW) power consumption at load nodes (rows: time steps, columns: load nodes); Nothing in forward-only mode
    power_producer::Union{Vector{Float64}, Nothing} # (MW) power output at producer node; Nothing in forward-only mode
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
    temperature::Union{Float64, Nothing}
end
```

# Fields
- `mass_flow`: total mass flow injected into the network in kg/s.
- `temperature`: producer outlet (supply) temperature in °C.
  `nothing` is only valid in `:backward_only` mode, where the forward thermal
  step is skipped and the producer supply temperature is not needed.

# Usage
The `policy` passed to [`run_simulation`](@ref) must return a `ProducerOutput`:

```julia
# 3-argument form (used in :full mode)
function policy(t, Tₐ, T_back)
    return ProducerOutput(mass_flow=15.0, temperature=90.0)
end

# 2-argument form (used in :forward_only, :backward_only, :hybrid modes)
function policy(t, Tₐ)
    return ProducerOutput(mass_flow=15.0, temperature=90.0)
end

# backward-only: temperature may be nothing
policy_bwd = [ProducerOutput(m_flow[i], nothing) for i in 1:N]
```
"""
struct ProducerOutput
    mass_flow::Float64
    temperature::Union{Float64, Nothing}   # Nothing only valid in :backward_only mode
end
ProducerOutput(; mass_flow, temperature=nothing) = ProducerOutput(mass_flow, temperature)

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
    field_value = getproperty(sr, s)
    if isnothing(field_value)
        return nothing
    end
    if occursin("load", String(s))
        return field_value[:, idx]
    else
        return field_value
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
run_simulation(network, sim_time, policy;
               mode=:full, T_return_inject=nothing,
               T0_f=60.0, T0_b=25.0, ambient_temperature=nothing)
```

This is the main entry point for time stepping.

REPEAT for N time steps:
1. Computes a steady-state hydraulic solution (mass flow distribution).
2. Advances thermal dynamics using the plug-flow method, following the
   steps enabled by `mode` (see below).

# Simulation Modes

| Step                | `:full` | `:forward_only` | `:backward_only` | `:hybrid` |
|---------------------|---------|-----------------|------------------|-----------|
| Hydraulics          | ✓       | ✓               | ✓                | ✓         |
| Forward thermal     | ✓       | ✓               | —                | ✓         |
| Load-model power    | ✓       | —               | —                | —         |
| Inject T_return     | —       | —               | ✓                | ✓         |
| Experimental power  | —       | —               | —                | ✓         |
| Backward thermal    | ✓       | —               | ✓                | ✓         |

- **`:full`** — complete simulation (default, unchanged behaviour).
- **`:forward_only`** — supply-pipe transport only; validates pipe delays and heat losses.
- **`:backward_only`** — return-pipe transport using injected measured return temperatures;
  validates return pipes independently of the load model.
- **`:hybrid`** — forward thermal + injected return temperatures + experimental power;
  decouples pipe transport from the load model.

# Output
- `SimulationResults` struct with time series of temperatures, flows, and powers.

# Arguments
- `network::Network`: prepared network (producer/load nodes identified, pipes attached).
- `sim_time`: equally spaced time vector.
    - `Vector{Float64}`: time in seconds.
    - `Vector{DateTime}`: timestamps (Δt is interpreted in seconds).
- `policy`: producer setpoints, either:
    - `Function` — called each step:
        - 3-arg `policy(t, Tₐ, T_back)` in `:full` mode.
        - 2-arg `policy(t, Tₐ)` in `:forward_only`, `:backward_only`, `:hybrid` modes.
    - `Vector{ProducerOutput}` — pre-built vector of length N.
      `temperature` may be `nothing` only in `:backward_only` mode.

# Keyword Arguments
- `mode`: simulation mode (`:full`, `:forward_only`, `:backward_only`, or `:hybrid`). Default `:full`.
- `T_return_inject`: `Dict{String, Vector{Float64}}` mapping load labels to injected return
  temperatures [°C] at each time step. Required for `:backward_only` and `:hybrid` modes.
- `T0_f`: initial temperature in the forward (supply) pipes [°C]. Default 60.0.
- `T0_b`: initial temperature in the backward (return) pipes [°C]. Default 25.0.
- `ambient_temperature`: optional `Vector{Float64}` of ambient temperatures [°C], length N.

# Returns
- `SimulationResults`: time series of temperatures, flows, and powers.

# Notes
- In `:hybrid` mode, `power_load` may contain negative values when the injected return
  temperature exceeds the simulated supply temperature at a load. This can occur due to
  timing mismatches, sensor noise, or calibration offsets. Negative values are preserved
  as diagnostics: a systematic negative bias signals a forward-pipe delay error.
- In `:backward_only` with `temperature=nothing`, `T_producer_out` is NaN and
  `power_producer` is `nothing`.
- The network structure is validated once at the start via `check_network!`.
- Time steps must be equally spaced.

# Chaining example
```julia
# Hybrid run: get simulated T_load_out with real return temperatures
sr_hyb = run_simulation(network, t, policy_vec; mode=:hybrid, T_return_inject=measured_T_v)

# Extract T_load_out as inject for a backward-only follow-up
T_inj2 = Dict(l => sr_hyb[l, :T_load_out] for l in keys(sr_hyb[:load_labels_dict]))
sr_bwd = run_simulation(network, t, policy_bwd; mode=:backward_only, T_return_inject=T_inj2)
```
"""
function run_simulation(
        network  :: Network,
        sim_time :: Union{Vector{Float64}, Vector{DateTime}},
        policy   :: Union{Function, Vector{ProducerOutput}};
        mode                :: Symbol = :full,
        T_return_inject     :: Union{Dict{String, Vector{Float64}}, Nothing} = nothing,
        T0_f                :: Float64 = 60.0,
        T0_b                :: Float64 = 25.0,
        ambient_temperature :: Union{Vector{Float64}, Nothing} = nothing) :: SimulationResults

    check_network!(network)

    # ---- parse time vector ----
    if eltype(sim_time) <: Float64
        dt = diff(sim_time)
    elseif eltype(sim_time) <: DateTime
        dt = diff(sim_time) .|> Dates.Second .|> Dates.value
    else
        error("Unsupported time vector element type: $(eltype(sim_time)). Expected Float64 or DateTime.")
    end
    if !all(isapprox.(dt, dt[1], atol=1e-10))
        error("Time steps are not equally spaced!")
    end
    Δt = float(dt[1])
    N  = length(sim_time)

    # ---- validate inputs ----
    mode ∈ (:full, :forward_only, :backward_only, :hybrid) ||
        error("Invalid mode: :$mode. Must be one of :full, :forward_only, :backward_only, :hybrid.")

    if mode ∈ (:backward_only, :hybrid)
        isnothing(T_return_inject) &&
            error("T_return_inject is required for mode=:$mode.")
        for (k, v) in T_return_inject
            k ∈ network.load_labels ||
                error("T_return_inject key \"$k\" is not a valid load label in the network.")
            length(v) == N ||
                error("T_return_inject[\"$k\"] has length $(length(v)), expected $N.")
        end
    end

    if policy isa Vector{ProducerOutput}
        length(policy) == N ||
            error("policy vector has length $(length(policy)), expected $N (length of sim_time).")
        if mode ∈ (:full, :forward_only, :hybrid)
            any(p -> isnothing(p.temperature), policy) &&
                error("policy Vector contains ProducerOutput with temperature=nothing, which is not allowed in mode=:$mode.")
        end
    else
        # test function policy for the first time step
        try
            Tₐ_test  = isnothing(ambient_temperature) ? nothing : ambient_temperature[1]
            test_out = if mode == :full
                policy(sim_time[1], Tₐ_test, T0_b)
            else
                policy(sim_time[1], Tₐ_test)
            end
            test_out isa ProducerOutput ||
                error("Policy function must return a ProducerOutput struct.")
            if mode ∈ (:full, :forward_only, :hybrid) && isnothing(test_out.temperature)
                error("Policy function returned ProducerOutput with temperature=nothing, which is not allowed in mode=:$mode.")
            end
        catch e
            error("Error calling policy for the first time step: ", e)
        end
    end

    # ---- initialise ----
    fill_pipes_with_initial_temperature!(network, T0_f, T0_b)

    num_loads        = length(network.load_labels)
    load_labels_cols = Dict(label => i for (i, label) in enumerate(network.load_labels))

    clamped_loads = Set{String}()  # tracks loads where return T was clamped (for a single end-of-run warning)

    results_mass_flow_producer       = Vector{Float64}(undef, N)
    results_mass_flow_load           = Matrix{Float64}(undef, N, num_loads)
    results_temperature_producer_out = fill(NaN, N)
    results_temperature_load_in      = fill(NaN, N, num_loads)
    results_temperature_producer_in  = mode == :forward_only ? nothing : fill(NaN, N)
    results_temperature_load_out     = mode == :forward_only ? nothing : fill(NaN, N, num_loads)
    results_power_consumption        = mode ∈ (:forward_only, :backward_only) ? nothing : fill(NaN, N, num_loads)

    # ---- time loop ----
    for i in 1:N
        Tₐ = isnothing(ambient_temperature) ? nothing : ambient_temperature[i]

        # policy dispatch
        input = if policy isa Vector{ProducerOutput}
            policy[i]
        elseif mode == :full
            T_back = i > 1 ? results_temperature_producer_in[i-1] : T0_b
            policy(sim_time[i], Tₐ, T_back)
        else
            policy(sim_time[i], Tₐ)
        end

        results_temperature_producer_out[i] = isnothing(input.temperature) ? NaN : input.temperature
        results_mass_flow_producer[i] = input.mass_flow

        # hydraulics (always)
        steady_state_hydronynamics!(network, input.mass_flow)

        # return_plugs collects cooled load plugs to feed into the backward step
        return_plugs = Dict{String, Plug}()

        # forward thermal (all modes except :backward_only)
        if mode != :backward_only
            output_plugs = time_step_thermal_dynamics_forward!(network, Δt, input.temperature, Tₐ)
            for (load_label, plug) in output_plugs
                col = load_labels_cols[load_label]
                results_temperature_load_in[i, col] = plug.T
                results_mass_flow_load[i, col]      = plug.m / Δt
            end
            return_plugs = output_plugs
        else
            # backward_only: read mass flows from network nodes; T_load_in is undefined (NaN)
            for (label, col) in load_labels_cols
                results_mass_flow_load[i, col] = network[label].common.mass_flow
            end
        end

        # load-model power (:full only)
        if mode == :full
            for load_label in keys(return_plugs)
                col = load_labels_cols[load_label]
                P = power_consumption(network[load_label], Tₐ)
                results_power_consumption[i, col] = P / 1000.0   # kW
                consume_power!(return_plugs[load_label], P, Δt) && push!(clamped_loads, load_label)
                results_temperature_load_out[i, col] = return_plugs[load_label].T
            end
        end

        # inject measured return temperatures (:backward_only and :hybrid)
        if mode ∈ (:backward_only, :hybrid)
            for (label, col) in load_labels_cols
                haskey(T_return_inject, label) || continue
                m_mass = results_mass_flow_load[i, col] * Δt
                T_inj  = T_return_inject[label][i]
                return_plugs[label] = Plug(T_inj, m_mass)
                results_temperature_load_out[i, col] = T_inj
                if mode == :hybrid
                    T_in  = results_temperature_load_in[i, col]
                    m_dot = results_mass_flow_load[i, col]
                    # P_exp [kW] = ṁ·cₚ·(T_load_in − T_load_out_inject)
                    # May be negative when the injected return temperature exceeds the simulated
                    # supply temperature at a load.  This can occur due to:
                    #   • timing mismatches between the simulation grid and the measurement grid
                    #     (resampling artefacts),
                    #   • measurement noise and sensor calibration offsets,
                    #   • real but short-lived thermal events not captured by the load model.
                    # Negative values are preserved (not clamped) because they are diagnostic:
                    # a systematic negative bias signals a forward-pipe delay error.
                    results_power_consumption[i, col] = m_dot * WATER_SPECIFIC_HEAT * (T_in - T_inj) / 1000.0
                end
            end
            # In :backward_only every load needs a return plug for the backward step.
            # Loads absent from T_return_inject (no measurement available) fall back to
            # T0_b so the backward thermal step always receives a complete plug dict.
            # Their T_load_out entries remain NaN and they are excluded from evaluation.
            if mode == :backward_only
                for (label, col) in load_labels_cols
                    haskey(return_plugs, label) && continue
                    return_plugs[label] = Plug(T0_b, results_mass_flow_load[i, col] * Δt)
                end
            end
        end

        # backward thermal (all modes except :forward_only)
        if mode != :forward_only && !isempty(return_plugs)
            incoming_plug = time_step_thermal_dynamics_backward!(network, Δt, return_plugs, Tₐ)
            results_temperature_producer_in[i] = incoming_plug.T
        end
    end

    # ---- warn once about clamped return temperatures ----
    if !isempty(clamped_loads)
        labels_str = join(sort(collect(clamped_loads)), ", ")
        @warn "Return temperature was clamped to the minimum ($(MINIMAL_RETURN_TEMPERATURE) °C) at one or more time steps for load(s): $labels_str"
    end

    # ---- producer power [MW] ----
    power_producer = if mode == :forward_only
        nothing
    elseif mode == :backward_only && all(isnan, results_temperature_producer_out)
        nothing
    else
        @. (results_temperature_producer_out[2:end] - results_temperature_producer_in[1:end-1]) *
           results_mass_flow_producer[1:end-1] * WATER_SPECIFIC_HEAT / 1_000_000.0
    end

    return SimulationResults(
        time               = sim_time,
        mass_flow_load     = results_mass_flow_load,
        mass_flow_producer = results_mass_flow_producer,
        T_load_in          = results_temperature_load_in,
        T_load_out         = results_temperature_load_out,
        T_producer_in      = results_temperature_producer_in,
        T_producer_out     = results_temperature_producer_out,
        power_load         = results_power_consumption,
        power_producer     = power_producer,
        load_labels_dict   = load_labels_cols
    )
end

# ------------------------------------------------- #
# Simulation helper methods
# ------------------------------------------------- #

"""Compute and assign relative mass-flow split coefficients `m_rel` on edges.

This performs a post-order traversal from leaves to root and sets, for each edge
leading into a node, the sum of `m_rel` values required downstream.

This is an internal step of [`steady_state_hydronynamics!`](@ref).
"""
function set_relative_mass_flows!(nw::Network)
    # iterate over nodes from leaves to root and set on each edge relative mass flow coefficient m_rel
    # we will do iterative post-order DFS traversal

    # first check that all leaf nodes have m_rel set
    for node in nw.load_labels
        if outdegree(nw, node) == 0 # leaf node
            if ismissing(nw[node].m_rel)
                error("Leaf node $(nw[node].common.info) has undefined m_rel. Please set m_rel on all load nodes before calling steady_state_hydrodynamics!")
            end
        end
    end

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

"""Assign absolute mass flows throughout the network.

Given the producer mass flow `mass_flow_source` [kg/s] and relative edge split
coefficients (see [`set_relative_mass_flows!`](@ref)), this propagates mass flows
from root to leaves and writes `mass_flow` to both nodes and pipe edges.
"""
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

"""Compute steady-state hydraulics for a `Network`.

This updates the mass-flow distribution throughout the network for a given
producer (source) mass flow `mass_flow_source` [kg/s].

Internally it:

1. computes relative flow splits (`set_relative_mass_flows!`),
2. assigns absolute mass flows on edges and nodes (`set_absolute_mass_flows!`).
"""
function steady_state_hydronynamics!(nw::Network, mass_flow_source::Float64)
    set_relative_mass_flows!(nw)
    set_absolute_mass_flows!(nw, mass_flow_source)
end

"""Initialize all pipe plug queues with uniform temperatures.

Fills every `InsulatedPipe` in the network with a single plug in the forward
queue at `temperature_f` and a single plug in the backward queue at
`temperature_b`.
"""
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
            if m_total > 0
                push!(edge.plugs_f, Plug(temperature_f, m_total))
                push!(edge.plugs_b, Plug(temperature_b, m_total))
            end
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

"""Collect plugs that exit a pipe over one time step.

Given a plug queue `plugs` (front = pipe outlet), mass flow `mass_flow` [kg/s],
and time step `Δt` [s], this pops and (if needed) splits plugs so that the
returned vector has total mass approximately `mass_flow*Δt`.
"""
function collect_exiting_water_plugs!(plugs::Vector{Plug}, mass_flow::Float64, Δt::Float64)::Vector{Plug}
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
function push_in_water_plugs_forward!(e::ZeroPipe, plugs::Vector{Plug})
    append!(e.plugs_f, plugs)
end
function push_in_water_plugs_backward!(e::ZeroPipe, plugs::Vector{Plug})
    append!(e.plugs_b, plugs)
end

"""Merge consecutive plugs whose temperatures are nearly equal.

This simplifies a plug queue in-place by combining adjacent plugs when
`abs(ΔT) < tol`.
"""
function merge_same_temperature_plugs!(plugs::Vector{Plug}; tol::Float64=1e-3)
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

"""Combine multiple plugs into a single mass-weighted average plug."""
function combine_plugs(plugs::Vector{Plug})::Plug
    @assert !isempty(plugs)
    total_mass = sum(p.m for p in plugs)
    if total_mass == 0.0
        return Plug(25.0, 0.0) # default temperature for zero mass
    end
    avg_temp = sum(p.T * p.m for p in plugs) / total_mass
    return Plug(avg_temp, total_mass)
end

"""Set the relative mass-flow coefficient for one load node."""
function set_load_params!(nw::Network, load_label::String, m_rel::Float64)
    @assert nw[load_label] isa LoadNode
    nw[load_label].m_rel = m_rel
end

"""Set the relative mass-flow coefficient for multiple load nodes.

`load_params` maps `load_label => m_rel`.
"""
function set_load_params!(nw::Network, load_params::Dict{String, Float64})
    for (load_label, m_rel) in load_params
        set_load_params!(nw, load_label, m_rel)
    end
end

"""Update the demand function parameters for a single load node, keeping the existing function.

Validates the new parameters with the node's current function before storing.
Raises an error if the node has no load function set yet — use [`set_load_fn!`](@ref) first.

See also: [`set_load_fn!`](@ref), [`validate_load_spec`](@ref).
"""
function set_load_params!(nw::Network, label::String, params::AbstractVector{<:Real};
                          T_a_range=-30.0:1.0:30.0)
    has_label(nw, label) || error("Node with label $label does not exist in the network.")
    node = nw[label]
    node isa LoadNode || error("Node with label $label is not a load node.")
    ismissing(node.load) && error("No load function set on node $label. Use set_load_fn! first.")
    p = convert(Vector{Float64}, params)
    validate_load_spec(node.load.fn, p; T_a_range=T_a_range)
    node.load.params = p
    nw[label] = node
end

"""Update demand function parameters for multiple load nodes, keeping each node's existing function.

`params_dict` maps `load_label => new_params`. All entries are validated before any node is updated.

See also: [`set_load_fn!`](@ref), [`set_load_params!`](@ref).
"""
function set_load_params!(nw::Network, params_dict::Dict{String, <:AbstractVector{<:Real}};
                          T_a_range=-30.0:1.0:30.0)
    for (label, params) in params_dict
        set_load_params!(nw, label, params; T_a_range=T_a_range)
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

"""Compute load power demand as a function of outdoor temperature.

Calls `node.load.fn(node.load.params, Tₐ)` and returns the result in **Watts**
(the load function returns kW — conversion is handled internally).
The result is clamped to zero to prevent negative power (energy flowing back into the network).
"""
function power_consumption(node::LoadNode, Tₐ::Float64)::Float64
    return max(0.0, node.load.fn(node.load.params, Tₐ)) * 1000.0
end
power_consumption(node::LoadNode, Tₐ::Nothing) = 0.0


"""Reduce a plug temperature by consuming `power` over a time step.

`power` is in Watts and `Δt` is in seconds. Updates `p` in-place.

Returns `true` if the return temperature was clamped to `MINIMAL_RETURN_TEMPERATURE`,
`false` otherwise.
"""
function consume_power!(p::Plug, power::Float64, Δt::Float64)::Bool
    cₚ = WATER_SPECIFIC_HEAT  # specific heat capacity in J/(kg·K)
    ΔT = power * Δt / (p.m * cₚ)  # temperature drop in K
    if p.T - ΔT < MINIMAL_RETURN_TEMPERATURE
        p.T = MINIMAL_RETURN_TEMPERATURE
        return true
    end
    p.T -= ΔT
    return false
end

"""Merge multiple plug sequences into one sequence (return-side merging).

This is used when multiple return streams meet at a junction: it produces a
single plug sequence that is consistent with the combined mass.
"""
function merge_water_plug_vectors!(plug_vectors::Vector{Vector{Plug}})::Vector{Plug}

    if length(plug_vectors) == 1
        merge_same_temperature_plugs!(plug_vectors[1]; tol=1e-2)
        return plug_vectors[1] # no need to merge if there is only one vector
    end

    @assert !any(isempty.(plug_vectors))


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

    # filter change points so there are not two too close to each other
    atol_kg = 1e-3 # tolerance of 1g ... that is the smallest plug vector that we will merge
    max_total_mass = maximum(total_masses)


    # now combine plugs between change points by mass-weighted average of temperature
    merged_plugs = Vector{Plug}()
    t_prev = 0.0
    for t in change_points

        # we dont want to work with this small plugs
        if (t-t_prev)*max_total_mass < atol_kg
            continue
        end

        if all(isempty.(plug_vectors)) # we already added the small drops to the last one interval
            break
        end

        plugs_at_t = Vector{Plug}()
        for i in eachindex(plug_vectors)
            target_mass = (t-t_prev) * total_masses[i]
            acc_mass = 0
            while acc_mass < target_mass && !isempty(plug_vectors[i])
                p = popfirst!(plug_vectors[i])
                
                if acc_mass + p.m > target_mass && p.m - (target_mass - acc_mass) >= atol_kg
                    # only part of the plug is in this interval
                    pushfirst!(plug_vectors[i], Plug(p.T, p.m - (target_mass - acc_mass))) # put the remaining part back to the front of the queue
                    p = Plug(p.T, (target_mass - acc_mass)) # take only the part of the plug that is in this interval
                end
                push!(plugs_at_t, p)
                acc_mass += p.m
            end
        end
        t_prev = t

        if any(isempty.(plug_vectors)) # if any of the branches in empty, add the remaining little drops if there are any
            # add all the other remaining small plugs at the end
            for i in eachindex(plug_vectors)
                while !isempty(plug_vectors[i])
                    push!(plugs_at_t, popfirst!(plug_vectors[i]))
                end
            end
            @assert all(isempty.(plug_vectors))
        end
        
        push!(merged_plugs, combine_plugs(plugs_at_t)) # combine plugs at this change point into one plug
    end

    # there should no plugs remain in the original vectors
    if !all(isempty, plug_vectors)
        @show plug_vectors
    end

    merge_same_temperature_plugs!(merged_plugs; tol=1e-2) # simplify the result by merging plugs with almost the same temperature
    
    @assert isapprox(sum(p.m for p in merged_plugs), sum(total_masses); atol=1e-6) # check that total mass is conserved
    return merged_plugs
end

"""Advance thermal dynamics by one time step.

This is a convenience wrapper used for manual stepping outside of [`run_simulation`](@ref).
It performs:

1. forward (supply) plug advection from producer to loads,
2. load power consumption and cooling,
3. backward (return) advection back to the producer.

# Arguments
- `nw::Network`: the network (must have steady-state mass flows already computed, e.g. via [`steady_state_hydronynamics!`](@ref)).
- `Δt::Float64`: time step in seconds.
- `input::ProducerOutput`: producer setpoints for this step.

# Keyword Arguments
- `ambient_temperature`: outdoor temperature in °C, or `nothing`. When `nothing`, load power consumption is skipped (loads don't cool the water) and pipe heat losses are not applied.

# Returns
- `(output_plugs, incoming_plug)` where `output_plugs` maps load labels to their inlet plug, and `incoming_plug` represents the return temperature entering the producer.
"""
function time_step_thermal_dynamics!(nw::Network, Δt::Float64, input::ProducerOutput; ambient_temperature::Union{Float64, Nothing}=nothing)
    output_plugs = time_step_thermal_dynamics_forward!(nw, Δt, input.temperature, ambient_temperature)

    if !isnothing(ambient_temperature)
        Tₐ_load = ambient_temperature
        clamped = String[]
        for load_label in keys(output_plugs)
            P = power_consumption(nw[load_label], Tₐ_load)
            consume_power!(output_plugs[load_label], P, Δt) && push!(clamped, load_label)
        end
        if !isempty(clamped)
            @warn "Return temperature clamped to minimum ($(MINIMAL_RETURN_TEMPERATURE) °C) at load(s): $(join(sort(clamped), ", "))"
        end
    end

    incoming_plug = time_step_thermal_dynamics_backward!(nw, Δt, output_plugs, ambient_temperature)
    return output_plugs, incoming_plug
end