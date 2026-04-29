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
        mass_flow_sump::Matrix{Float64}
        T_sump_f::Matrix{Float64}
        T_sump_b::Union{Matrix{Float64}, Nothing}
        sump_labels::Dict{String, Int}
        producer_label::String
        m_rel_load::Matrix{Float64}
        tau1_load::Matrix{Float64}
        tau2_load::Matrix{Float64}
        tau1_producer::Vector{Float64}
        tau2_producer::Vector{Float64}
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
- `mass_flow_sump`: sump mass flows in kg/s. Size `N × nsumps`.
- `T_sump_f`: supply (forward) temperature at each sump in °C. Size `N × nsumps`. `NaN` in backward-only mode.
- `T_sump_b`: return (backward) temperature at each sump in °C. Size `N × nsumps`. `Nothing` in forward-only mode.
- `sump_labels`: mapping from sump label to column index used in the `*_sump` matrices.
- `producer_label`: label of the producer node (used by the universal indexing aliases).
- `m_rel_load`: relative mass-flow split coefficients at load nodes. Size `N × nloads`. `m_rel_load[i, j]` is the scalar `m_rel` value used at load `j` during time step `i`. For constant loads this column is uniform; for time-varying loads it captures the full trajectory.
- `tau1_load`: front-wall (oldest exiting fluid) supply-path transit delay at each load [s]. Size `N × nloads`. `NaN` for steps where initial-fill water is still exiting (before the pipe is fully flushed).
- `tau2_load`: back-wall (newest exiting fluid) supply-path transit delay at each load [s]. Size `N × nloads`. `NaN` for initial-fill steps.
- `tau1_producer`: front-wall return-path transit delay at the producer [s]. Length `N`. `NaN` for initial-fill steps or in `:forward_only` mode.
- `tau2_producer`: back-wall return-path transit delay at the producer [s]. Length `N`. `NaN` for initial-fill steps or in `:forward_only` mode.

# Indexing
Convenience accessors are provided:

- `sr[:time]` returns the time vector.
- `sr[:load_labels]` returns the load labels.
- `sr[:load_labels_dict]` returns the label→column dictionary.
- `sr[:sump_labels]` returns the sump labels.
- `sr[:sump_labels_dict]` returns the sump label→column dictionary.
- `sr["L1", :T_load_in]` returns the time series for load L1 (a vector).
- `sr["L1", :m_rel_load]` returns the m_rel time series for load L1 (a vector).
- `sr["S1", :T_sump_f]` returns the supply temperature time series for sump S1 (a vector).
- `sr["S1", :T_sump_b]` returns the return temperature time series for sump S1 (a vector).
- `sr["S1", :mass_flow_sump]` returns the mass flow time series for sump S1 (a vector).
- `sr["L1", :tau1_load]` returns the front-wall delay time series for load L1 (a vector, NaN = initial fill).
- `sr["L1", :tau2_load]` returns the back-wall delay time series for load L1 (a vector, NaN = initial fill).
- `sr[:tau1_producer]` returns the front-wall return-path delay vector at the producer.
- `sr[:tau2_producer]` returns the back-wall return-path delay vector at the producer.

Universal aliases (work for any label — load, sump, or producer):

- `sr["L1", :T_in]` returns the inlet temperature: supply T for loads/sumps, return T for the producer.
- `sr["L1", :T_out]` returns the outlet temperature: return T for loads/sumps, supply T for the producer.
- `sr["L1", :mass_flow]` returns the mass flow time series for the named node.

# Notes
- All matrices are organized as `(time step, node index)`.
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
    mass_flow_sump::Matrix{Float64}         # mass flows at sump nodes (rows: time steps, columns: sump nodes)
    T_sump_f::Matrix{Float64}              # supply (forward) temperature at sump nodes; NaN in backward_only mode
    T_sump_b::Union{Matrix{Float64}, Nothing}  # return (backward) temperature at sump nodes; Nothing in forward_only mode
    sump_labels::Dict{String, Int}          # labels of sump nodes corresponding to columns
    producer_label::String                  # label of the producer node
    m_rel_load::Matrix{Float64}             # relative mass-flow coefficients at load nodes (rows: time steps, columns: load nodes)
    tau1_load::Matrix{Float64}              # front-wall supply-path delay at each load [s]; NaN = initial fill
    tau2_load::Matrix{Float64}              # back-wall supply-path delay at each load [s]; NaN = initial fill
    tau1_producer::Vector{Float64}          # front-wall return-path delay at producer [s]; NaN = initial fill or forward_only
    tau2_producer::Vector{Float64}          # back-wall return-path delay at producer [s]; NaN = initial fill or forward_only
end
# constructor with keyword arguments for better readability
function SimulationResults(;time, mass_flow_load, mass_flow_producer, T_load_in, T_load_out, T_producer_in, T_producer_out, power_load, power_producer, load_labels_dict, mass_flow_sump, T_sump_f, T_sump_b, sump_labels_dict, producer_label, m_rel_load, tau1_load, tau2_load, tau1_producer, tau2_producer)
    return SimulationResults(time, mass_flow_load, mass_flow_producer, T_load_in, T_load_out, T_producer_in, T_producer_out, power_load, power_producer, load_labels_dict, mass_flow_sump, T_sump_f, T_sump_b, sump_labels_dict, producer_label, m_rel_load, tau1_load, tau2_load, tau1_producer, tau2_producer)
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

# Constructor
```julia
ProducerOutput(; mass_flow, temperature=nothing)
```

# Fields
- `mass_flow`: total mass flow rate injected into the network **in kg/s** (a continuous
  rate, independent of the time-step size Δt). The hydraulic solver uses this value
  directly to distribute flow across all pipes.
- `temperature`: producer outlet (supply) temperature in °C.
  `nothing` is only valid in `:backward_only` mode, where the forward thermal
  step is skipped and the producer supply temperature is not needed.

# Usage
The `policy` passed to [`run_simulation`](@ref) must return a `ProducerOutput`:

```julia
# Unified 3-argument form — always called as policy(t, T_a, T_back).
# T_a   is missing when no ambient_temperature vector is provided.
# T_back is missing in :forward_only, :backward_only, and :hybrid modes.
function policy(t, T_a, T_back)
    T_back_eff = ismissing(T_back) ? 40.0 : T_back   # fallback when return T unavailable
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
    # Sump-specific fields
    sump_syms = (:T_sump_f, :T_sump_b, :mass_flow_sump)
    if s in sump_syms
        haskey(sr.sump_labels, label) || error("Sump label \"$label\" not found in SimulationResults.")
        idx = sr.sump_labels[label]
        field_value = getproperty(sr, s)
        isnothing(field_value) && return nothing
        return field_value[:, idx]
    end

    # Load fields
    options = setdiff(fieldnames(SimulationResults), [:load_labels, :sump_labels, :mass_flow_sump, :T_sump_f, :T_sump_b])
    if(s ∉ options && s != :load_labels_dict)
        error("$(s): Invalid symbol for SimulationResults getindex. Valid symbols are: $(options), or sump symbols: $(sump_syms)")
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
    options = (fieldnames(SimulationResults)..., :load_labels_dict, :sump_labels_dict)
    if(s ∉ options)
        error("$(s): Invalid symbol for SimulationResults getindex. Valid symbols are: $(options)")
    end
    if s == :load_labels
        return collect(keys(sr.load_labels))
    elseif s == :sump_labels
        return collect(keys(sr.sump_labels))
    elseif s == :load_labels_dict
        return sr.load_labels
    elseif s == :sump_labels_dict
        return sr.sump_labels
    end
    return getproperty(sr, s)
end


Base.length(sr::SimulationResults) = length(sr.time)
Base.show(io::IO, sr::SimulationResults) = print(io, "SimulationResults with $(length(sr)) time steps, $(length(sr.load_labels)) load node(s) and $(length(sr.sump_labels)) sump node(s).")

"""
k₀ is the first index of output at load with `label`, that doesn't contain any initial-fill water.
"""
function get_k₀(sr::SimulationResults, mode::Symbol=:fwd; label::Union{Nothing, String, Vector{String}}=nothing)::Union{Nothing, Int, Dict{String, Union{Nothing, Int}}}
    if mode ∉ (:fwd, :bwd)
        ArgumentError("Invalid mode: $mode. Must be :fwd or :bwd.")
    end
    if mode == :fwd
        if isnothing(label)
            label = sr[:load_labels] # all load labels
        end
        if label isa String
            if label ∉ sr[:load_labels]
                error("Label \"$label\" not found in SimulationResults load labels.")
            else
                label = [label] # wrap in vector for uniform processing
            end
        end
        # label is now a Vector{String}
        k₀_dict = Dict{String, Union{Nothing, Int}}()
        for l in label
            if mode == :fwd
                # this is the first backwall, that comes from the producer.
                # Because frontwall is also the backwall of the previous plug, this exiting plug is a mixture of initial fill and producer water.
                tau2_mixture = findfirst(!isnan, sr[l, :tau2_load])
                k₀_dict[l] = isnothing(tau2_mixture) ? nothing : tau2_mixture + 1 # we have to take +1, because k₀ doesn't contain any initial-fill water
            else # mode == :bwd
                # this is the first frontwall, that comes from the load
                tau1_mixture = findfirst(!isnan, sr[l, :tau1_load])
                k₀_dict[l] = isnothing(tau1_mixture) ? nothing : tau1_mixture + 1 # we have to take +1, because k₀ doesn't contain any initial-fill water
            end
        end
        if length(label) == 1
            return first(values(k₀_dict)) # return single Int if only one label was requested
        else
            return k₀_dict # return Dict if multiple labels were requested
        end
    else # mode == :bwd, we have only one vector of delays, label is ignored
        if !isnothing(label)
            @warn "Label argument is ignored in :bwd mode, because there is only one producer return temperature time series."
        end
        tau1_mixture = findfirst(!isnan, sr[:tau1_producer])
        return isnothing(tau1_mixture) ? nothing : tau1_mixture + 1 # we have to take +1, because k₀ doesn't contain any initial-fill water
    end
end
get_k0 = get_k₀  # alias with ASCII character


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
- `policy`: producer setpoints for each time step, either:
    - `Function` — always called as `policy(t, T_a, T_back)`, must return a
      [`ProducerOutput`](@ref):
        - `T_a` is `missing` when no `ambient_temperature` vector is supplied.
        - `T_back` is `missing` in `:forward_only`, `:backward_only`, and `:hybrid` modes;
          it holds the previous-step producer return temperature in `:full` mode.
    - `Vector{ProducerOutput}` — pre-built vector of length N.
      `temperature` may be `nothing` only in `:backward_only` mode.
    In both cases, `ProducerOutput.mass_flow` must be in **kg/s** (a continuous flow rate,
    not kg per time step).

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

    # ---- validate m_rel mode uniformity and vector lengths ----
    has_const_m_rel  = any(l -> !(network[l].m_rel isa Vector{Float64}), network.load_labels)
    has_vector_m_rel = any(l ->   network[l].m_rel isa Vector{Float64},  network.load_labels)
    if has_const_m_rel && has_vector_m_rel
        error("m_rel is constant on some load nodes and time-varying on others. " *
              "All loads must use the same form (constant Float64 or Vector{Float64}) uniformly.")
    end
    if has_vector_m_rel
        for label in network.load_labels
            m = network[label].m_rel
            if m isa Vector{Float64} && length(m) != N
                error("Load \"$label\": m_rel vector has length $(length(m)), " *
                      "but sim_time has length $N. They must match.")
            end
        end
    end

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
            Tₐ_test     = isnothing(ambient_temperature) ? missing : ambient_temperature[1]
            T_back_test = (mode == :full) ? T0_b : missing
            test_out    = policy(sim_time[1], Tₐ_test, T_back_test)
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

    # Pre-compute pipe m_rel coefficients once.
    # For constant loads this writes a Float64 scalar to each pipe (single DFS pass).
    # For time-varying loads this writes a Vector{Float64} of length N to each pipe (N DFS passes).
    # Inside the loop only set_absolute_mass_flows! is called, reading m_rel(pipe, i).
    set_relative_mass_flows!(network)

    num_loads        = length(network.load_labels)
    load_labels_cols = Dict(label => i for (i, label) in enumerate(network.load_labels))

    num_sumps        = length(network.sump_labels)
    sump_labels_cols = Dict(label => i for (i, label) in enumerate(network.sump_labels))

    clamped_loads = Set{String}()  # tracks loads where return T was clamped (for a single end-of-run warning)

    results_mass_flow_producer       = Vector{Float64}(undef, N)
    results_mass_flow_load           = Matrix{Float64}(undef, N, num_loads)
    results_temperature_producer_out = fill(NaN, N)
    results_temperature_load_in      = fill(NaN, N, num_loads)
    results_temperature_producer_in  = mode == :forward_only ? nothing : fill(NaN, N)
    results_temperature_load_out     = mode == :forward_only ? nothing : fill(NaN, N, num_loads)
    results_power_consumption        = mode ∈ (:forward_only, :backward_only) ? nothing : fill(NaN, N, num_loads)
    results_mass_flow_sump           = Matrix{Float64}(undef, N, num_sumps)
    results_T_sump_f                 = fill(NaN, N, num_sumps)
    results_T_sump_b                 = mode == :forward_only ? nothing : fill(NaN, N, num_sumps)
    results_m_rel_load               = Matrix{Float64}(undef, N, num_loads)
    results_tau1_load                = fill(NaN, N, num_loads)
    results_tau2_load                = fill(NaN, N, num_loads)
    results_tau1_producer            = fill(NaN, N)
    results_tau2_producer            = fill(NaN, N)

    # ---- time loop ----
    for i in 1:N
        Tₐ = isnothing(ambient_temperature) ? nothing : ambient_temperature[i]

        # policy dispatch
        input = if policy isa Vector{ProducerOutput}
            policy[i]
        else
            Tₐ_policy     = isnothing(Tₐ) ? missing : Tₐ
            T_back_policy = (mode == :full) ? (i > 1 ? results_temperature_producer_in[i-1] : T0_b) : missing
            policy(sim_time[i], Tₐ_policy, T_back_policy)
        end

        results_temperature_producer_out[i] = isnothing(input.temperature) ? NaN : input.temperature
        results_mass_flow_producer[i] = input.mass_flow

        # hydraulics: m_rel already precomputed; only propagate absolute flows
        set_absolute_mass_flows!(network, input.mass_flow, i)

        # record sump mass flows and load m_rel values (always, after hydraulics)
        for (label, col) in sump_labels_cols
            results_mass_flow_sump[i, col] = network[label].common.mass_flow
        end
        for (label, col) in load_labels_cols
            results_m_rel_load[i, col] = _load_m_rel(network[label].m_rel, i)
        end

        # return_plugs collects cooled load plugs to feed into the backward step
        return_plugs = Dict{String, Plug}()

        # forward thermal (all modes except :backward_only)
        if mode != :backward_only
            sump_plugs_f = Dict{String, Plug}()
            output_plugs = time_step_thermal_dynamics_forward!(network, Δt, i, input.temperature, Tₐ; sump_plugs=sump_plugs_f)
            for (load_label, plug) in output_plugs
                col = load_labels_cols[load_label]
                results_temperature_load_in[i, col] = plug.T
                results_mass_flow_load[i, col]      = plug.m / Δt
                isnan(plug.t1) || (results_tau1_load[i, col] = (i - 1 - plug.t1) * Δt)
                isnan(plug.t2) || (results_tau2_load[i, col] = (i     - plug.t2) * Δt)
            end
            for (label, col) in sump_labels_cols
                haskey(sump_plugs_f, label) && (results_T_sump_f[i, col] = sump_plugs_f[label].T)
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
                P = power_consumption(network[load_label], Tₐ, i)
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
            sump_plugs_b = !isnothing(results_T_sump_b) ? Dict{String, Plug}() : nothing
            incoming_plug = time_step_thermal_dynamics_backward!(network, Δt, i, return_plugs, Tₐ; sump_plugs=sump_plugs_b)
            results_temperature_producer_in[i] = incoming_plug.T
            isnan(incoming_plug.t1) || (results_tau1_producer[i] = (i - 1 - incoming_plug.t1) * Δt)
            isnan(incoming_plug.t2) || (results_tau2_producer[i] = (i     - incoming_plug.t2) * Δt)
            if !isnothing(sump_plugs_b)
                for (label, col) in sump_labels_cols
                    haskey(sump_plugs_b, label) && (results_T_sump_b[i, col] = sump_plugs_b[label].T)
                end
            end
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
        load_labels_dict   = load_labels_cols,
        mass_flow_sump     = results_mass_flow_sump,
        T_sump_f           = results_T_sump_f,
        T_sump_b           = results_T_sump_b,
        sump_labels_dict   = sump_labels_cols,
        producer_label     = network.producer_label::String,
        m_rel_load         = results_m_rel_load,
        tau1_load          = results_tau1_load,
        tau2_load          = results_tau2_load,
        tau1_producer      = results_tau1_producer,
        tau2_producer      = results_tau2_producer,
    )
end

# ------------------------------------------------- #
# Simulation helper methods
# ------------------------------------------------- #

# Extract the scalar m_rel value for the current time step from a LoadNode's m_rel field.
# Constant Float64 → returned as-is; Vector{Float64} → indexed by step.
_load_m_rel(m_rel::Float64, ::Int) = m_rel
_load_m_rel(m_rel::Vector{Float64}, step::Int) = m_rel[step]

"""Compute and assign relative mass-flow split coefficients `m_rel` on edges.

Two overloads are provided:

---

    set_relative_mass_flows!(nw, step::Int)

**Per-step** form. Performs one post-order DFS traversal and writes a `Float64`
scalar to the `m_rel` field of every pipe edge. `step` selects which element of
a time-varying `LoadNode.m_rel` vector to use; for constant `Float64` loads it
is ignored.

Use this inside a manual time-stepping loop when you prefer to recompute the
split at every step rather than pre-allocating vectors.

---

    set_relative_mass_flows!(nw)

**Vectorised / precompute** form. Called once *before* a simulation loop.

- If all load nodes carry a constant `Float64` `m_rel`: performs a single DFS
  and writes `Float64` scalars to the pipe edges (same as `step=1`).
- If all load nodes carry a time-varying `Vector{Float64}` `m_rel` of length N:
  pre-allocates a `Vector{Float64}(undef, N)` on every pipe, then runs N DFS
  passes (one per step) filling `pipe.m_rel[k]` at each pass.

After calling this form, use [`m_rel(pipe, step)`](@ref) to read the scalar for
any given step without re-running the DFS.

---

Both overloads **normalise** load `m_rel` values so they sum to `1.0` across all
load nodes at each time step before propagating upstream.  Trunk-pipe `m_rel`
values are therefore proper fractions: a junction pipe carries the sum of its
subtree's normalised leaf fractions, and the producer's outgoing pipe always
carries `m_rel = 1.0`.

Both forms are internal steps of [`steady_state_hydrodynamics!`](@ref) and are
called automatically by [`run_simulation`](@ref).

See also: [`set_absolute_mass_flows!`](@ref), [`steady_state_hydrodynamics!`](@ref).
"""
function set_relative_mass_flows!(nw::Network, step::Union{Int, Nothing}=nothing)
    for label in nw.load_labels
        if outdegree(nw, label) == 0 && ismissing(nw[label].m_rel)
            error("Leaf node $(nw[label].common.info) has undefined m_rel. Please set m_rel on all load nodes before calling set_relative_mass_flows!")
        end
    end

    has_vector = any(l -> nw[l].m_rel isa Vector{Float64}, nw.load_labels)

    if !isnothing(step)
        _dfs_m_rel!(nw, step)
        return
    end

    if !has_vector
        _dfs_m_rel!(nw, 1)
        return
    end

    N = maximum(length(nw[l].m_rel) for l in nw.load_labels if nw[l].m_rel isa Vector{Float64})

    # Per-step normalisation totals: totals[k] = sum of all leaf load m_rel at step k
    totals = zeros(N)
    for l in nw.load_labels
        outdegree(nw, l) == 0 && (totals .+= nw[l].m_rel)
    end

    # Pre-allocate Vector{Float64}(undef, N) on every pipe
    for (src_lbl, dst_lbl) in MetaGraphsNext.edge_labels(nw.mg)
        pipe = nw[src_lbl, dst_lbl]
        if pipe isa InsulatedPipe || pipe isa ZeroPipe
            pipe.m_rel = Vector{Float64}(undef, N)
        end
    end

    # Single post-order DFS: fill entire m_rel vectors in one pass
    root = nw.producer_label::String
    stack = Vector{Tuple{String, Bool}}()
    push!(stack, (root, false))
    while !isempty(stack)
        node, visited = pop!(stack)
        if visited
            nw[node] isa ProducerNode && continue
            parent_node = inneighbors(nw, node)[1]
            pipe = nw[parent_node, node]
            if outdegree(nw, node) == 0
                @assert nw[node] isa LoadNode "Leaf node $(nw[node].common.info) is not a LoadNode. Please check the network structure."
                pipe.m_rel .= nw[node].m_rel ./ totals
            elseif nw[node] isa JunctionNode || nw[node] isa SumpNode
                fill!(pipe.m_rel, 0.0)
                for child in outneighbors(nw, node)
                    pipe.m_rel .+= nw[node, child].m_rel
                end
            end
        else
            push!(stack, (node, true))
            for child in outneighbors(nw, node)
                push!(stack, (child, false))
            end
        end
    end
end

# Write m_rel to a pipe: scalar assignment or vector-element write depending on current type.
function _pipe_m_rel_write!(pipe::Union{InsulatedPipe, ZeroPipe}, step::Int, val::Float64)
    if pipe.m_rel isa Vector{Float64}
        pipe.m_rel[step] = val
    else
        pipe.m_rel = val
    end
end

# One post-order DFS pass that assigns m_rel on every pipe edge for `step`.
# Load m_rel values are normalised to sum to 1.0 across all load nodes before propagating,
# so trunk pipes carry proper fractions and the root pipe always carries 1.0.
function _dfs_m_rel!(nw::Network, step::Int)
    root = nw.producer_label::String
    total_m_rel = sum(_load_m_rel(nw[l].m_rel, step) for l in nw.load_labels if outdegree(nw, l) == 0)
    stack = Vector{Tuple{String, Bool}}()
    push!(stack, (root, false))
    while !isempty(stack)
        node, visited = pop!(stack)
        if visited
            nw[node] isa ProducerNode && continue
            parent_node = inneighbors(nw, node)[1]
            pipe = nw[parent_node, node]
            if outdegree(nw, node) == 0
                @assert nw[node] isa LoadNode "Leaf node $(nw[node].common.info) is not a LoadNode. Please check the network structure."
                _pipe_m_rel_write!(pipe, step, _load_m_rel(nw[node].m_rel, step) / total_m_rel)
            elseif nw[node] isa JunctionNode || nw[node] isa SumpNode
                m_sum = sum(m_rel(nw[node, child], step) for child in outneighbors(nw, node))::Float64
                _pipe_m_rel_write!(pipe, step, m_sum)
            end
        else
            push!(stack, (node, true))
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

The optional `step` argument (default `1`) is passed to [`m_rel(pipe, step)`](@ref)
to read the split coefficient for this time step. It is relevant when pipes carry
pre-computed time-varying `Vector{Float64}` `m_rel` (written by the no-argument
[`set_relative_mass_flows!`](@ref) overload). For scalar `Float64` pipe `m_rel`,
`step` is ignored.
"""
function set_absolute_mass_flows!(nw::Network, mass_flow_source::Float64, step::Int=1)
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

        sum_rel_flow = sum(m_rel(nw[node, child], step) for child in outneighbors(nw, node))
        if sum_rel_flow == 0.0
            warn("Node $(nw[node].common.info) has zero total relative mass flow to its children. Skipping mass flow assignment for its outgoing edges.")
            continue # no flow through this node
        end
        for child in outneighbors(nw, node)
            edge = nw[node, child]
            edge.mass_flow = nw[node].common.mass_flow * m_rel(edge, step) / sum_rel_flow
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

The optional `step` argument (default `1`) is forwarded to [`set_relative_mass_flows!`](@ref)
and selects which element of each time-varying `m_rel` vector to use (see [`LoadNode`](@ref)).
For constant `m_rel` loads the value is ignored. When called outside of [`run_simulation`](@ref)
(e.g. for a single hydraulic snapshot), the default `step=1` is always correct for constant loads.
"""
function steady_state_hydrodynamics!(nw::Network, mass_flow_source::Float64, step::Int=1)
    set_relative_mass_flows!(nw, step)
    set_absolute_mass_flows!(nw, mass_flow_source, step)
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

function time_step_thermal_dynamics_forward!(nw::Network, Δt::Float64, step::Int, temperature_source::Float64, ambient_temperature::Union{Float64, Nothing}=nothing; sump_plugs::Union{Dict{String,Plug}, Nothing}=nothing)::Dict{String, Plug}
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
        new_plug = Plug(temperature_source, plug_mass, Float64(step) - 0.5, Float64(step - 1), Float64(step))
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

        # collect plugs exiting from parent edge; heat loss is applied per-plug inside collect
        plugs = if !isnothing(ambient_temperature) && parent_edge isa InsulatedPipe
            collect_exiting_water_plugs!(parent_edge.plugs_f, parent_edge.mass_flow, Δt, step;
                d=inner_diameter(parent_edge), R=heat_resistance_forward(parent_edge), T_a=ambient_temperature)
        else
            collect_exiting_water_plugs!(parent_edge.plugs_f, parent_edge.mass_flow, Δt, step)
        end

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

        # record supply temperature at sump nodes
        if !isnothing(sump_plugs) && nw[node] isa SumpNode
            sump_plugs[node] = combine_plugs(plugs)
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
                    next_plug = Plug(p.T, next_mass, p.k, p.t1, p.t2) # inherit k, t1, t2
                    push!(next_pipe_plugs, next_plug)
                end
            end
            # push the smaller plugs into the child edge
            push_in_water_plugs_forward!(edge, next_pipe_plugs)
        end
    end

    return output_plugs
end

function time_step_thermal_dynamics_backward!(nw::Network, Δt::Float64, step::Int, incoming_plugs::Dict{String, Plug}, ambient_temperature::Union{Float64, Nothing}=nothing; sump_plugs::Union{Dict{String,Plug}, Nothing}=nothing)::Union{Plug, Nothing}
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
                plug = incoming_plugs[node]
                plug.k  = Float64(step) - 0.5  # midpoint of step: return plug injected uniformly during this step
                plug.t1 = Float64(step - 1)
                plug.t2 = Float64(step)
                push_in_water_plugs_backward!(parent_edge, [plug])

            else # visiting for second time, so all children have been processed
                children_plug_vectors = Vector{Vector{Plug}}() # collect plugs from child edges to merge them
                # collect plugs exiting from child edges and apply heat loss
                for child in outneighbors(nw, node)
                    edge = nw[node, child]
                    # collect plugs; heat loss applied per-plug inside collect
                    plugs = if !isnothing(ambient_temperature) && edge isa InsulatedPipe
                        collect_exiting_water_plugs!(edge.plugs_b, edge.mass_flow, Δt, step;
                            d=inner_diameter(edge), R=heat_resistance_backward(edge), T_a=ambient_temperature)
                    else
                        collect_exiting_water_plugs!(edge.plugs_b, edge.mass_flow, Δt, step)
                    end
                    push!(children_plug_vectors, plugs)
                end
                # combine plugs from all child edges; k_avg assigned from interval midpoints
                merged = merge_water_plug_vectors!(children_plug_vectors, step)
                # record return temperature at sump nodes
                if !isnothing(sump_plugs) && nw[node] isa SumpNode
                    sump_plugs[node] = combine_plugs(merged)
                end
                # push merged plugs back to parent edge
                if node != root
                    parent = inneighbors(nw, node)[1] # assume there is only one parent
                    parent_edge = nw[parent, node]
                    push_in_water_plugs_backward!(parent_edge, merged) # k already set by merge_water_plug_vectors!
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

    return output_plug
end

"""Collect plugs that exit a pipe over one time step.

Given a plug queue `plugs` (front = pipe outlet), mass flow `mass_flow` [kg/s],
time step `Δt` [s], and current simulation step `step`, this pops and (if needed)
splits plugs so that the returned vector has total mass approximately `mass_flow*Δt`.

Each returned plug has its `k` field set to the fractional step at which the plug's
mass midpoint exits the pipe — this is the entry time into the next pipe segment.
The transit time through the current pipe is `(k_new - k_old) * Δt`, giving
sub-step resolution.

If `d`, `R`, and `T_a` are provided, `apply_exit_heat_loss!` is called per plug
before its `k` is updated, using the continuous transit time `(k_exit - k_entry)·Δt`.
"""
function collect_exiting_water_plugs!(
    plugs::Vector{Plug}, mass_flow::Float64, Δt::Float64, step::Int;
    d::Union{Float64, Nothing}=nothing,
    R::Union{Float64, Nothing}=nothing,
    T_a::Union{Float64, Nothing}=nothing
)::Vector{Plug}
    exited_plugs = Vector{Plug}()
    M_exit = mass_flow * Δt
    if M_exit <= 0.0
        return exited_plugs
    end
    apply_heat_loss = !isnothing(d) && !isnothing(R) && !isnothing(T_a)
    M_acc = 0.0
    while !isempty(plugs) && M_acc < M_exit
        p = popfirst!(plugs)
        if M_acc + p.m > M_exit
            # only part of the plug exits; split into exiting and remaining portions
            m_exits   = M_exit - M_acc
            m_remains = p.m - m_exits
            k_entry    = p.k
            k_avg_exit = (step - 1) + (M_acc + m_exits / 2) / M_exit
            # split [t1, t2] proportionally by mass: front portion exits, back portion remains
            frac    = m_exits / p.m
            t_split = p.t1 + (p.t2 - p.t1) * frac
            # remaining (back) portion stays in the pipe with the original entry k
            pushfirst!(plugs, Plug(p.T, m_remains, k_entry, t_split, p.t2))
            # exiting (front) portion: apply heat loss, then update k
            exiting = Plug(p.T, m_exits, k_avg_exit, p.t1, t_split)
            apply_heat_loss && apply_exit_heat_loss!(exiting, k_entry, k_avg_exit, d, R, Δt, T_a)
            push!(exited_plugs, exiting)
            M_acc = M_exit
            break
        else
            # entire plug exits; t1/t2 unchanged, only k is updated
            k_entry    = p.k
            k_avg_exit = (step - 1) + (M_acc + p.m / 2) / M_exit
            apply_heat_loss && apply_exit_heat_loss!(p, k_entry, k_avg_exit, d, R, Δt, T_a)
            p.k = k_avg_exit
            push!(exited_plugs, p)
            M_acc += p.m
        end
    end
    if isempty(plugs) && !isapprox(M_acc, M_exit; atol=1e-6)
        error("Not enough mass in the pipe to exit! Mass exited: $M_acc, expected: $M_exit. This indicates a problem with the simulation, possibly due to numerical errors or incorrect mass flow assignment.")
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

"""Combine multiple plugs into a single mass-weighted average plug.

Temperature and k are mass-weighted averages. t1 is taken from the first plug
(oldest/front boundary) and t2 from the last plug (newest/back boundary),
preserving the full injection time window of the combined interval.
"""
function combine_plugs(plugs::Vector{Plug})::Plug
    @assert !isempty(plugs)
    total_mass = sum(p.m for p in plugs)
    if total_mass == 0.0
        return Plug(25.0, 0.0, 0.0, NaN, NaN)
    end
    avg_temp = sum(p.T * p.m for p in plugs) / total_mass
    avg_k    = sum(p.k * p.m for p in plugs) / total_mass
    return Plug(avg_temp, total_mass, avg_k, plugs[1].t1, plugs[end].t2)
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
    validate_load_spec(node.load.fn, p; T_a_range=T_a_range,
                       use_mass_flow=node.load.use_mass_flow, use_time=node.load.use_time)
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

# Apply heat loss to a single plug exiting an InsulatedPipe.
# τ = (k_exit - k_entry) * Δt is the continuous transit time through the pipe.
# τ_c = ρ·cₚ·A·R is the thermal time constant of the pipe [s].
function apply_exit_heat_loss!(p::Plug, k_entry::Float64, k_exit::Float64, d::Float64, R::Float64, Δt::Float64, T_a::Float64)
    τ = (k_exit - k_entry) * Δt
    if τ > 0.0
        A   = 1/4 * π * d^2
        τ_c = WATER_DENSITY * WATER_SPECIFIC_HEAT * A * R # [s]
        p.T = T_a + (p.T - T_a) * exp(-τ / τ_c)
    end
end

"""Compute load power demand for a single time step.

Dispatches to the load function with the arguments dictated by `node.load.use_mass_flow`
and `node.load.use_time` (see [`LoadSpec`](@ref)).

- `Tₐ::Float64` — ambient temperature in °C.
- `step::Int` — 1-based simulation step index, passed to the function when `use_time = true`.

Returns power in **Watts** (the load function returns kW — conversion is handled internally).
The result is clamped to zero to prevent negative power (energy flowing back into the network).
"""
function power_consumption(node::LoadNode, Tₐ::Float64, step::Int=1)::Float64
    v = if node.load.use_mass_flow && node.load.use_time
        m = mass_flow(node)
        node.load.fn(node.load.params, Tₐ, m, step)
    elseif node.load.use_mass_flow
        m = mass_flow(node)
        node.load.fn(node.load.params, Tₐ, m)
    elseif node.load.use_time
        node.load.fn(node.load.params, Tₐ, step)
    else
        node.load.fn(node.load.params, Tₐ)
    end
    if ismissing(v)
        @warn "Load function for node $(node.common.info) returned missing value at step $step. Treating as zero power consumption."
        return 0.0   # out-of-range step (e.g. lookup vector shorter than sim length)
    end
    return max(0.0, Float64(v)) * 1000.0
end
power_consumption(node::LoadNode, Tₐ::Nothing, ::Int=1) = 0.0


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
function merge_water_plug_vectors!(plug_vectors::Vector{Vector{Plug}}, step::Int)::Vector{Plug}

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
    atol_kg = 1e-3 # tolerance of 1g ... that is the smallest plug we will merge
    max_total_mass = maximum(total_masses)

    # now combine plugs between change points by mass-weighted average of temperature;
    # k_avg for each merged plug is the midpoint of its interval on the common f ∈ [0,1] axis.
    merged_plugs = Vector{Plug}()
    t_prev = 0.0
    for t in change_points

        # skip negligible intervals
        if (t - t_prev) * max_total_mass < atol_kg
            continue
        end

        if all(isempty.(plug_vectors)) # we already added the small drops to the last interval
            break
        end

        plugs_at_t = Vector{Plug}()
        for i in eachindex(plug_vectors)
            target_mass = (t - t_prev) * total_masses[i]
            acc_mass = 0.0
            while acc_mass < target_mass && !isempty(plug_vectors[i])
                p = popfirst!(plug_vectors[i])
                if acc_mass + p.m > target_mass && p.m - (target_mass - acc_mass) >= atol_kg
                    # only part of this plug falls in the current interval; split t1/t2 proportionally
                    portion = target_mass - acc_mass
                    frac    = portion / p.m
                    t_split = p.t1 + (p.t2 - p.t1) * frac
                    pushfirst!(plug_vectors[i], Plug(p.T, p.m - portion, p.k, t_split, p.t2))
                    p = Plug(p.T, portion, p.k, p.t1, t_split)
                end
                push!(plugs_at_t, p)
                acc_mass += p.m
            end
        end

        t_interval_start = t_prev
        t_prev = t

        if any(isempty.(plug_vectors)) # if any branch is exhausted, absorb remaining tiny drops
            for i in eachindex(plug_vectors)
                while !isempty(plug_vectors[i])
                    push!(plugs_at_t, popfirst!(plug_vectors[i]))
                end
            end
            @assert all(isempty.(plug_vectors))
        end

        # mass-weighted temperature, k, t1, t2; k_avg = midpoint of interval in the step
        m_total_interval = sum(p.m for p in plugs_at_t)
        T_merged     = sum(p.T  * p.m for p in plugs_at_t) / m_total_interval
        k_avg_merged = Float64(step - 1) + (t_interval_start + t) / 2
        t1_merged    = sum(p.t1 * p.m for p in plugs_at_t) / m_total_interval
        t2_merged    = sum(p.t2 * p.m for p in plugs_at_t) / m_total_interval
        push!(merged_plugs, Plug(T_merged, m_total_interval, k_avg_merged, t1_merged, t2_merged))
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
- `nw::Network`: the network (must have steady-state mass flows already computed, e.g. via [`steady_state_hydrodynamics!`](@ref)).
- `Δt::Float64`: time step in seconds.
- `input::ProducerOutput`: producer setpoints for this step.

# Keyword Arguments
- `ambient_temperature`: outdoor temperature in °C, or `nothing`. When `nothing`, load power consumption is skipped (loads don't cool the water) and pipe heat losses are not applied.
- `step`: simulation step counter. Each plug's fractional entry time and exit time are derived from `step`, giving continuous transit times `τ = (k_exit - k_entry) · Δt`. Defaults to `1`. When calling this function in a manual stepping loop, pass the iteration index so that heat loss is computed correctly.

# Returns
- `(output_plugs, incoming_plug)` where `output_plugs` maps load labels to their inlet plug, and `incoming_plug` represents the return temperature entering the producer.
"""
function time_step_thermal_dynamics!(nw::Network, Δt::Float64, input::ProducerOutput; ambient_temperature::Union{Float64, Nothing}=nothing, step::Int=1)
    output_plugs = time_step_thermal_dynamics_forward!(nw, Δt, step, input.temperature, ambient_temperature)

    if !isnothing(ambient_temperature)
        Tₐ_load = ambient_temperature
        clamped = String[]
        for load_label in keys(output_plugs)
            P = power_consumption(nw[load_label], Tₐ_load, step)
            consume_power!(output_plugs[load_label], P, Δt) && push!(clamped, load_label)
        end
        if !isempty(clamped)
            @warn "Return temperature clamped to minimum ($(MINIMAL_RETURN_TEMPERATURE) °C) at load(s) $(length(clamped)): $(join(sort(clamped), ", "))"
        end
    end

    incoming_plug = time_step_thermal_dynamics_backward!(nw, Δt, step, output_plugs, ambient_temperature)
    return output_plugs, incoming_plug
end