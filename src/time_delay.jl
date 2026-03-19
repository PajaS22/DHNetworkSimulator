"""
    compute_time_delay(network, sim_result, Δt; label=nothing) -> τ_vec or Dict

Compute the transit time delay τ(t) for every simulation time step t.

τ(t) is defined as the time [s] it takes for a plug of water that **enters the
producer node at step t** to travel through the supply network and reach the
specified load node.  The plug traverses each pipe segment sequentially; its
travel time through a pipe equals the pipe's water volume divided by the (time-
varying) mass flow through that pipe.

# Arguments
- `network::Network`            — network whose topology is used (must have been
  set up with `identify_producer_and_loads!` and `fill_physical_params!`).
- `sim_result::SimulationResults` — completed simulation result; load mass flows
  are read from `sim_result[lbl, :mass_flow_load]` and used as pipe flows.
- `Δt::Real`                    — simulation time step [s].

# Keyword arguments
- `label::Union{String,Nothing}` — if a string, compute τ only for that load and
  return a single `Vector{Union{Float64,Missing}}`.  If `nothing` (default),
  compute for **all** loads and return a `Dict{String, Vector{Union{Float64,Missing}}}`.

# Output
The returned vector (or each vector in the dict) has the same length N as
`sim_result[:mass_flow_producer]`.  Entry τ[t] contains the delay in seconds,
or `missing` when the plug cannot reach the load within the simulation window
(i.e. insufficient mass-flow data remain after step t to fill the pipe path).

# Notes
- The flow through an intermediate pipe segment is computed as the **sum of the
  simulated mass flows** of all downstream loads (those reachable through that pipe).
- Only `InsulatedPipe` segments contribute to the delay; `ZeroPipe` and other
  edge types are skipped.
- The calculation uses a precomputed cumulative-sum array per pipe segment and
  `searchsortedfirst` for O(log N) binary search per time step.

# Example
```julia
sr = run_simulation(network, t_vec, policy; mode=:forward_only,
                    ambient_temperature=T_a_vec)
τ = compute_time_delay(network, sr, Δt; label="M2_VS_3")
τ_all = compute_time_delay(network, sr, Δt)   # Dict for every load
```
"""
function compute_time_delay(
    network::Network,
    sim_result::SimulationResults,
    Δt::Real;
    label::Union{String, Nothing} = nothing,
)::Union{Vector{Union{Float64, Missing}}, Dict{String, Vector{Union{Float64, Missing}}}}

    labels = label === nothing ? collect(network.load_labels) : [label]
    N      = length(sim_result[:mass_flow_producer])

    results = Dict{String, Vector{Union{Float64, Missing}}}()

    for lbl in labels
        path = _td_path_to_load(network, lbl)   # [producer, …, lbl]
        segs = _td_pipe_segments(network, sim_result, path, Float64(Δt), N)

        τ = Vector{Union{Float64, Missing}}(missing, N)

        for t in 1:N
            start   = t      # next pipe begins accumulating from this step
            success = true

            for (M_pipe, cumflow) in segs
                # cumflow[i] = Σ_{j=1}^{i-1} flow[j]*Δt  (length N+1)
                # We need: cumflow[k+1] - cumflow[start] >= M_pipe
                #       => cumflow[k+1] >= cumflow[start] + M_pipe
                start > N && (success = false; break)
                target = cumflow[start] + M_pipe
                # Binary search in cumflow[start+1 .. N+1]
                idx = searchsortedfirst(cumflow, target,
                                        start + 1, N + 1, Base.Order.Forward)
                idx > N + 1 && (success = false; break)
                start = idx - 1   # k = idx-1 is the step where pipe is flushed
            end

            success && (τ[t] = (start - t) * Δt)
        end

        results[lbl] = τ
    end

    return label === nothing ? results : results[label]
end


"""
    compute_initial_delay(network, sim_result, Δt; label=nothing)

Compute only τ(1) — the transit time [s] for a water plug that enters the
producer node at the **first** simulation step to reach each load.

This is the **warmup duration**: temperature signals recorded before this time
still reflect the water that pre-filled the pipes at simulation start, not the
actual producer output, and should therefore be excluded when comparing
simulated results to measurements.

Equivalent to `compute_time_delay(...)[1]` for each load, but runs in
O(P·log N) per load instead of O(N·P·log N), making it fast to call on large
networks before any analysis loop.

# Arguments
- `network::Network`: the network (topology + pipe geometry).
- `sim_result::SimulationResults`: completed simulation; load mass flows are
  read from `sim_result[lbl, :mass_flow_load]`.
- `Δt::Real`: simulation time step [s].

# Keyword arguments
- `label::Union{String,Nothing}`: if a string, compute only for that load and
  return a single `Union{Float64,Missing}`.  If `nothing` (default), compute
  for all loads and return a `Dict{String,Union{Float64,Missing}}`.

# Returns
`missing` when the plug does not reach the load within the simulation window.

# Example
```julia
sr = run_simulation(network, t_vec, policy; mode=:forward_only)
warmup = compute_initial_delay(network, sr, Δt)   # Dict for all loads

# Skip warmup steps when comparing to measurements at stride `s`:
skip(lbl) = ceil(Int, coalesce(warmup[lbl], 0.0) / Δt / s)
```
"""
function compute_initial_delay(
    network::Network,
    sim_result::SimulationResults,
    Δt::Real;
    label::Union{String, Nothing} = nothing,
)::Union{Union{Float64, Missing}, Dict{String, Union{Float64, Missing}}}

    labels = label === nothing ? collect(network.load_labels) : [label]
    N      = length(sim_result[:mass_flow_producer])

    results = Dict{String, Union{Float64, Missing}}()

    for lbl in labels
        path = _td_path_to_load(network, lbl)
        segs = _td_pipe_segments(network, sim_result, path, Float64(Δt), N)

        start   = 1   # plug enters at step 1
        success = true
        for (M_pipe, cumflow) in segs
            start > N && (success = false; break)
            target = cumflow[start] + M_pipe
            idx = searchsortedfirst(cumflow, target,
                                    start + 1, N + 1, Base.Order.Forward)
            idx > N + 1 && (success = false; break)
            start = idx - 1   # arrival step for this pipe
        end

        results[lbl] = success ? (start - 1) * Float64(Δt) : missing
    end

    return label === nothing ? results : results[label]
end


# ── Internal helpers ──────────────────────────────────────────────────────────

# Walk backwards from `load_label` to the producer and return the full node path
# [producer, …, load_label].  Assumes a tree (each node has at most one parent).
function _td_path_to_load(network::Network, load_label::String)::Vector{String}
    path = String[load_label]
    node = load_label
    while !isempty(inneighbors(network, node))
        parent = inneighbors(network, node)[1]
        pushfirst!(path, parent)
        node = parent
    end
    return path
end

# Return the set of load labels reachable (downstream) from `node`, including
# `node` itself if it is a load.
function _td_downstream_loads(network::Network, node::String)::Set{String}
    node in network.load_labels && return Set{String}([node])
    loads = Set{String}()
    for child in outneighbors(network, node)
        union!(loads, _td_downstream_loads(network, child))
    end
    return loads
end

# For each InsulatedPipe segment on `path`, build a tuple (M_pipe_kg, cumflow)
# where cumflow is the length-(N+1) precomputed cumulative integral of the
# mass flow [kg] through that pipe over time.
# Flow through an intermediate pipe = sum of simulated mass flows of all
# downstream loads (those reachable through the pipe's destination node).
function _td_pipe_segments(
    network::Network,
    sim_result::SimulationResults,
    path::Vector{String},
    Δt::Float64,
    N::Int,
)::Vector{Tuple{Float64, Vector{Float64}}}

    segs = Tuple{Float64, Vector{Float64}}[]
    sim_load_labels = sim_result[:load_labels]

    for i in 1:length(path) - 1
        src, dst = path[i], path[i + 1]
        pipe = network[src, dst]
        pipe isa InsulatedPipe || continue

        # Pipe water volume and mass
        r = inner_diameter(pipe) / 2.0
        V = π * r^2 * pipe_length(pipe)           # m³
        M = WATER_DENSITY * V                      # kg

        # Aggregate flow of all downstream loads through this pipe
        flow = zeros(Float64, N)
        for dl in _td_downstream_loads(network, dst)
            dl in sim_load_labels || continue
            flow .+= sim_result[dl, :mass_flow_load]
        end

        # Cumulative mass [kg] starting from time index 1
        # cumflow[i] = Σ_{j=1}^{i-1} flow[j]*Δt  =>  length N+1
        cumflow = Vector{Float64}(undef, N + 1)
        cumflow[1] = 0.0
        for j in 1:N
            cumflow[j + 1] = cumflow[j] + flow[j] * Δt
        end

        push!(segs, (M, cumflow))
    end

    return segs
end
