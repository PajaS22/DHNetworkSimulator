"""
    compute_time_delay(network, sim_result, Δt; label=nothing) -> τ_vec or Dict

Compute the transit time delay τ(t) for every simulation time step t.

τ(t) is defined as the time [s] the plug of water **exiting at step t** has
spent travelling through the supply network from the producer to the specified
load node.  Equivalently, `T_supply(t − τ(t))` is the producer temperature of
the fluid arriving at the load at step t.

The plug traverses each pipe segment sequentially; its travel time through a
pipe equals the pipe's water volume divided by the (time-varying) mass flow
through that pipe.

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
or `missing` when the pipe path was not yet fully flushed by step t (i.e.
insufficient mass-flow data precede step t to account for the pipe volume).

# Notes
- The flow through an intermediate pipe segment is computed as the **sum of the
  simulated mass flows** of all downstream loads (those reachable through that pipe).
- Only `InsulatedPipe` segments contribute to the delay; `ZeroPipe` and other
  edge types are skipped.
- The calculation uses a precomputed cumulative-sum array per pipe segment and
  `searchsortedfirst` for O(log N) binary search per time step.  Pipe segments
  are traversed in **reverse** (load → producer) to find the injection time of
  the plug exiting at each step.

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
            end_step = t      # plug exits the last pipe at step t
            success  = true

            for (M_pipe, cumflow) in Iterators.reverse(segs)
                # cumflow[i] = Σ_{j=1}^{i-1} flow[j]*Δt  (length N+1)
                # We need cumflow[s] <= cumflow[end_step+1] - M_pipe < cumflow[s+1]
                # => entry step s = searchsortedfirst(...) - 1
                target = cumflow[end_step + 1] - M_pipe
                target < 0 && (success = false; break)
                # Binary search in cumflow[1 .. end_step+1]
                idx = searchsortedfirst(cumflow, target,
                                        1, end_step + 1, Base.Order.Forward)
                end_step = idx - 1   # step at which plug entered this pipe
            end

            success && (τ[t] = (t - end_step) * Δt)
        end

        results[lbl] = τ
    end

    return label === nothing ? results : results[label]
end


"""
    compute_initial_delay(network, sim_result, Δt; label=nothing)

Compute the **warmup duration** [s] for each load — the time until the pipe
path is first fully flushed with producer-side fluid.

This equals τ(t₀), the transit delay at the first valid output step t₀ (the
earliest step at which a complete column of producer fluid has passed through
all pipe segments to reach the load).  Temperature signals before this time
still reflect pre-filled water and should be excluded when comparing simulated
results to measurements.

Equivalent to `compute_time_delay(...)[t₀]` for each load, but runs in
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
