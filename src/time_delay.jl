"""
    compute_time_delay(network, sim_result, Δt; label=nothing, delay_type=:back)

Compute the plug-flow transit delay for every simulation time step.

At each output step `t`, the exiting fluid was injected over a window of input
steps.  That window has two boundaries:

- **Front wall** — the oldest boundary (fluid that has been in the pipe the
  longest).  It entered the first pipe at a continuous time *within* some input
  step, and exits at the **start** of output step `t` (time `(t−1)·Δt`).
- **Back wall** — the newest boundary (fluid that just entered the exit window).
  It entered the first pipe within a later input step and exits at the **end**
  of output step `t` (time `t·Δt`).

`delay_type` selects which wall's transit time to return:

| `delay_type` | returned τ | formula (single pipe, Planning notation) |
|:---|:---|:---|
| `:front` | time front wall spent in pipe | `(dτ₁ − 1 + α) · Δt` |
| `:back`  | time back wall spent in pipe  | `(dτ₂ + 1 − β) · Δt` |
| `:both`  | named tuple `(; front, back)` | both vectors |

Both delays are measured as physical transit time (continuous entry → continuous
exit), so they are directly comparable with `τ₁` / `τ₂` from the Planning
module's `continuous_time_delays`.

# Arguments
- `network::Network`              — topology (must have been set up with
  `identify_producer_and_loads!` and `fill_physical_params!`).
- `sim_result::SimulationResults` — completed simulation; load mass flows are
  read from `sim_result[lbl, :mass_flow_load]`.
- `Δt::Real`                      — simulation time step [s].

# Keyword arguments
- `label::Union{String,Nothing}` — if a string, compute for that load only and
  return a vector (or named tuple for `:both`).  If `nothing` (default),
  compute for all loads and return a `Dict`.
- `delay_type::Symbol`           — `:front`, `:back` (default), or `:both`.

# Output
- `:front` / `:back` with a label → `Vector{Union{Float64,Missing}}`, length N.
- `:both` with a label → `NamedTuple` `(; front::Vector{…}, back::Vector{…})`.
- Without a label → `Dict{String, <above>}`.

Entry is `missing` when the pipe path was not yet fully flushed at step `t`.

# Notes
- Only `InsulatedPipe` segments contribute to the delay.
- Pipe segments are traversed in **reverse** (load → producer).  The fractional
  correction (α / β) is applied only at the outermost (producer-side) pipe,
  where the fluid originally entered the network.

# Example
```julia
sr = run_simulation(network, t_vec, policy; mode=:forward_only)
τ_back  = compute_time_delay(network, sr, Δt; label="L1")
τ_front = compute_time_delay(network, sr, Δt; label="L1", delay_type=:front)
τ_both  = compute_time_delay(network, sr, Δt; label="L1", delay_type=:both)
τ_both.front   # Vector{Union{Float64,Missing}}
τ_both.back    # Vector{Union{Float64,Missing}}
```
"""
function compute_time_delay(
    network::Network,
    sim_result::SimulationResults,
    Δt::Real;
    label::Union{String, Nothing} = nothing,
    delay_type::Symbol = :back,
)
    delay_type ∈ (:front, :back, :both) ||
        error("delay_type must be :front, :back, or :both, got :$delay_type")

    labels = label === nothing ? collect(network.load_labels) : [label]
    N      = length(sim_result[:mass_flow_producer])

    results = Dict{String, Any}()

    for lbl in labels
        path = path_to_load(network, lbl)
        segs = path_pipe_segments(network, sim_result, path, Float64(Δt), N)

        τ_front = (delay_type == :back)  ? nothing :
                  Vector{Union{Float64, Missing}}(missing, N)
        τ_back  = (delay_type == :front) ? nothing :
                  Vector{Union{Float64, Missing}}(missing, N)

        for t in 1:N
            # ── back wall: entered at continuous time within step end_step_b ──
            # target_b = M_new = cumflow[t+1] - M_pipe
            # The back wall is the newest fluid exiting at step t.
            # It entered the first pipe at fraction β through step end_step_b,
            # i.e. at time (end_step_b - 1 + β)·Δt, and exits at time t·Δt.
            # τ_back = (t - end_step_b + 1 - β)·Δt
            if !isnothing(τ_back)
                end_step = t
                success  = true
                last_cumflow = segs[end][2]   # will be overwritten; just a placeholder
                last_target  = 0.0

                for (M_pipe, cumflow) in Iterators.reverse(segs)
                    target = cumflow[end_step + 1] - M_pipe
                    target < 0 && (success = false; break)
                    idx      = searchsortedfirst(cumflow, target,
                                                 1, end_step + 1, Base.Order.Forward)
                    end_step     = idx - 1
                    last_cumflow = cumflow
                    last_target  = target
                end

                if success && end_step >= 1
                    denom = last_cumflow[end_step + 1] - last_cumflow[end_step]
                    β = denom > 0 ? (last_target - last_cumflow[end_step]) / denom : 0.0
                    τ_back[t] = (t - end_step + 1 - β) * Δt
                end
            end

            # ── front wall: entered at continuous time within step end_step_f ──
            # target_f = M_prev = cumflow[t] - M_pipe
            # The front wall is the oldest fluid still exiting at step t.
            # It entered at fraction (1-α) through step end_step_f,
            # i.e. at time (end_step_f - α)·Δt, and exits at time (t-1)·Δt.
            # τ_front = (t - end_step_f - 1 + α)·Δt
            if !isnothing(τ_front)
                end_step = t
                success  = true
                last_cumflow = segs[end][2]
                last_target  = 0.0

                for (M_pipe, cumflow) in Iterators.reverse(segs)
                    target = cumflow[end_step] - M_pipe
                    target < 0 && (success = false; break)
                    idx      = searchsortedfirst(cumflow, target,
                                                 1, end_step + 1, Base.Order.Forward)
                    end_step     = idx - 1
                    last_cumflow = cumflow
                    last_target  = target
                end

                if success && end_step >= 1
                    denom = last_cumflow[end_step + 1] - last_cumflow[end_step]
                    α = denom > 0 ? (last_cumflow[end_step + 1] - last_target) / denom : 1.0
                    τ_front[t] = (t - end_step - 1 + α) * Δt
                end
            end
        end

        results[lbl] = if delay_type == :both
            (; front=τ_front, back=τ_back)
        elseif delay_type == :front
            τ_front
        else
            τ_back
        end
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
        path = path_to_load(network, lbl)
        segs = path_pipe_segments(network, sim_result, path, Float64(Δt), N)

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


# ── Path and segment helpers ───────────────────────────────────────────────────

"""
    path_to_load(network, load_label) -> Vector{String}

Return the ordered node-label path `[producer, intermediate..., load_label]`
from the producer to `load_label`, walking backwards via `inneighbors`.
Assumes a tree topology (each node has at most one parent).
"""
function path_to_load(network::Network, load_label::String)::Vector{String}
    path = String[load_label]
    node = load_label
    while !isempty(inneighbors(network, node))
        parent = inneighbors(network, node)[1]
        pushfirst!(path, parent)
        node = parent
    end
    return path
end

"""
    downstream_loads(network, node) -> Set{String}

Return the set of load labels reachable downstream from `node`,
including `node` itself if it is a load node.
"""
function downstream_loads(network::Network, node::String)::Set{String}
    node in network.load_labels && return Set{String}([node])
    loads = Set{String}()
    for child in outneighbors(network, node)
        union!(loads, downstream_loads(network, child))
    end
    return loads
end

"""
    path_pipe_segments(network, sim_result, path, Δt, N)

For each `InsulatedPipe` on `path` (ordered node labels from producer to load),
return a `(M_pipe_kg, cumflow)` tuple where `M_pipe_kg` is the fluid mass [kg]
and `cumflow` is the length-`(N+1)` cumulative mass through the pipe, computed
from the simulated mass flows of all downstream loads.
"""
function path_pipe_segments(
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
        for dl in downstream_loads(network, dst)
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
