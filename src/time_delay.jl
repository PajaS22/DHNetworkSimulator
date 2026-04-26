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
For `:front`, non-missing at step `t` means the **oldest** fluid exiting at `t`
is traceable to producer input through every pipe — i.e., the entire exit
window is from producer input and the output is "clean".  For `:back`,
non-missing only means the *newest* fluid is traceable, which can happen while
the front wall still traces to pre-fill in an upstream pipe; use `:front` (or
[`compute_k0`](@ref) / [`compute_initial_delay`](@ref)) to identify the first
fully-flushed step.

# Notes
- Only `InsulatedPipe` segments contribute to the delay.
- Delays are computed **per pipe** using exact integer-step indices `dτ₁`/`dτ₂`
  and fractional corrections `α`/`β` (ported from the Planning module's
  `time_delays`), then **composed** by step-shifting from load back to producer.
  The fractional correction is taken from the **producer-side pipe** evaluated
  at the shifted entry step, giving continuous-time accuracy at every pipe boundary.

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
    dt     = Float64(Δt)

    results = Dict{String, Any}()

    for lbl in labels
        path        = path_to_load(network, lbl)
        pipe_delays = _path_pipe_delays(network, sim_result, path, dt, N)

        τ_front = (delay_type == :back)  ? nothing :
                  Vector{Union{Float64, Missing}}(missing, N)
        τ_back  = (delay_type == :front) ? nothing :
                  Vector{Union{Float64, Missing}}(missing, N)

        if isempty(pipe_delays)
            # No InsulatedPipe segments on path — zero delay everywhere.
            τ_back  === nothing || fill!(τ_back,  0.0)
            τ_front === nothing || fill!(τ_front, 0.0)
        else
            # pipe_delays_rev: load-side first, producer-side last.
            # Iterating load→producer lets us step-shift j backward through pipes.
            # The fractional correction (α/β) is read from the producer-side pipe
            # (pipe_delays_rev[end]) at the step j BEFORE its delay is applied.
            pipe_delays_rev = reverse(pipe_delays)

            for t in 1:N
                # ── back wall (newest plug, dτ₂, β) ────────────────────────────
                # τ₂[t] = (total_dτ₂ + 1 − β) · Δt
                # β comes from the producer-side pipe at its output step j₁,
                # which is the step j right before we apply p₁'s dτ₂ shift.
                if !isnothing(τ_back)
                    j = t; total_dτ = 0; ok = true; j_first = t
                    for pd in pipe_delays_rev
                        ismissing(pd.dτ₂[j]) && (ok = false; break)
                        j_first = j          # j before this pipe's shift
                        d       = Int(pd.dτ₂[j])
                        total_dτ += d
                        j -= d
                    end
                    if ok && !ismissing(pipe_delays_rev[end].βₖ[j_first])
                        β = Float64(pipe_delays_rev[end].βₖ[j_first])
                        τ_back[t] = (total_dτ + 1 - β) * dt
                    end
                end

                # ── front wall (oldest plug, dτ₁, α) ───────────────────────────
                # τ₁[t] = (total_dτ₁ − 1 + α) · Δt
                # α comes from the producer-side pipe at its output step j₁.
                if !isnothing(τ_front)
                    j = t; total_dτ = 0; ok = true; j_first = t
                    for pd in pipe_delays_rev
                        ismissing(pd.dτ₁[j]) && (ok = false; break)
                        j_first = j
                        d       = Int(pd.dτ₁[j])
                        total_dτ += d
                        j -= d
                    end
                    if ok && !ismissing(pipe_delays_rev[end].αₖ[j_first])
                        α = Float64(pipe_delays_rev[end].αₖ[j_first])
                        τ_front[t] = (total_dτ - 1 + α) * dt
                    end
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

Compute the **warmup duration** [s] for each load — the transit time of the
front wall (oldest exiting fluid) at the first step where the pipe path output
is entirely from producer-side fluid.

This equals `τ_front[k₀]`, where `k₀` is the first output step at which the
**front wall** (oldest fluid in the exit window) can be traced back to
producer input through every pipe on the path.  Because the front wall is the
oldest exiting fluid, a non-missing front-wall delay at step `k₀` guarantees
that *all* fluid exiting at `k₀` is from producer input — the back wall is
necessarily from producer input too, since it entered later.  Relying on the
back-wall delay alone (`τ_back[k₀_back]`) gives an earlier, potentially wrong
`k₀` because the back wall (newest fluid) can be from producer input while the
front wall (oldest fluid) still traces back to pre-fill in an upstream pipe.

Temperature signals before step `k₀` still reflect pre-filled water and should
be excluded when comparing simulated results to measurements.

Computed with full fractional precision using the same per-pipe `dτ₁`/`dτ₂`/`α`/`β`
method as [`compute_time_delay`](@ref).

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

    τ_front_all = compute_time_delay(network, sim_result, Δt;
                                     label=label, delay_type=:front)

    if label !== nothing
        k₀ = findfirst(!ismissing, τ_front_all)
        return k₀ === nothing ? missing : τ_front_all[k₀]
    else
        results = Dict{String, Union{Float64, Missing}}()
        for (lbl, τ_vec) in τ_front_all
            k₀ = findfirst(!ismissing, τ_vec)
            results[lbl] = k₀ === nothing ? missing : τ_vec[k₀]
        end
        return results
    end
end


"""
    compute_k0(network, sim_result, Δt; label=nothing)

Return the **first meaningful output step index** `k₀` for each load — the
earliest step at which the entire exit window is from producer-side fluid.

`k₀` is the first step at which the **front wall** (oldest exiting fluid) can
be traced back to producer input through every pipe on the path.  Because the
front wall is the oldest fluid in the exit window, a non-missing front-wall
delay at `k₀` guarantees that *all* fluid exiting at `k₀` is from producer
input.  Steps `1 : k₀-1` contain pre-fill fluid and should be excluded from
analysis.

Using the back wall instead of the front wall would give an earlier, potentially
wrong `k₀`: the back wall (newest fluid) can arrive from producer input while
the front wall (oldest fluid) still traces back to pre-fill in an upstream pipe.

# Arguments
- `network::Network`: network topology and pipe geometry.
- `sim_result::SimulationResults`: completed simulation.
- `Δt::Real`: simulation time step [s].

# Keyword arguments
- `label::Union{String,Nothing}`: if a string, compute only for that load and
  return a single `Union{Int,Missing}`.  If `nothing` (default), compute for
  all loads and return a `Dict{String,Union{Int,Missing}}`.

# Returns
`missing` when the pipe path is never fully flushed within the simulation window.

# Example
```julia
sr  = run_simulation(network, t_vec, policy; mode=:forward_only)
k₀  = compute_k0(network, sr, Δt; label="L1")   # Int — first valid step
τ₀  = compute_initial_delay(network, sr, Δt; label="L1")   # Float64 [s]

# Analyse only valid steps:
T_valid = sr["L1", :T_load_in][k₀:end]

# Equivalent one-liner for all loads:
k0_dict = compute_k0(network, sr, Δt)
```

See also: [`compute_initial_delay`](@ref), [`compute_time_delay`](@ref).
"""
function compute_k0(
    network::Network,
    sim_result::SimulationResults,
    Δt::Real;
    label::Union{String, Nothing} = nothing,
)::Union{Union{Int, Missing}, Dict{String, Union{Int, Missing}}}

    τ_front_all = compute_time_delay(network, sim_result, Δt;
                                     label=label, delay_type=:front)

    if label !== nothing
        k₀ = findfirst(!ismissing, τ_front_all)
        return k₀ === nothing ? missing : k₀
    else
        results = Dict{String, Union{Int, Missing}}()
        for (lbl, τ_vec) in τ_front_all
            k₀ = findfirst(!ismissing, τ_vec)
            results[lbl] = k₀ === nothing ? missing : k₀
        end
        return results
    end
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


# ── Private per-pipe delay computation ────────────────────────────────────────

# Compute plug-flow time delays for a single pipe.
# Returns a named tuple (dτ₁, dτ₂, αₖ, βₖ, k₀) of length-N vectors.
# Ported from the Planning module's time_delays() in pipe_transport.jl.
#
# Conventions (1-based, matching Planning):
#   cum_m[k] = Σᵢ₌₁ᵏ ṁ[i]·Δt  — cumulative mass through end of step k
#   k₀ = first output step where all exiting fluid entered during steps 1..N
#   dτ₁[k] = integer delay to oldest plug (front wall) exiting at step k
#   dτ₂[k] = integer delay to newest plug (back wall) exiting at step k
#   αₖ[k]  = fraction of the dτ₁ input plug exiting at step k (for τ₁ formula)
#   βₖ[k]  = fraction of the dτ₂ input plug exiting at step k (for τ₂ formula)
#
# Continuous-time delays follow:
#   τ₁[k] = (dτ₁[k] - 1 + αₖ[k]) · Δt   (front-wall transit time)
#   τ₂[k] = (dτ₂[k] + 1 - βₖ[k]) · Δt   (back-wall transit time)
function _pipe_time_delays(
    Δt::Float64,
    M_pipe::Float64,
    ṁ::Vector{Float64},
)
    N     = length(ṁ)
    cum_m = cumsum(ṁ) .* Δt   # cum_m[k] = total mass entering pipe through end of step k

    # First output step where all exiting fluid comes from the simulation window [1, N]
    k₀_idx = findfirst(x -> x > M_pipe, cum_m)
    if k₀_idx === nothing
        # Pipe never fully flushed within N steps
        return (
            dτ₁ = Vector{Union{Int,Missing}}(missing, N),
            dτ₂ = Vector{Union{Int,Missing}}(missing, N),
            αₖ  = Vector{Union{Float64,Missing}}(missing, N),
            βₖ  = Vector{Union{Float64,Missing}}(missing, N),
            k₀  = N + 1,
        )
    end

    dτ₁ = Vector{Union{Int,Missing}}(missing, N)
    dτ₂ = Vector{Union{Int,Missing}}(missing, N)
    αₖ  = Vector{Union{Float64,Missing}}(missing, N)
    βₖ  = Vector{Union{Float64,Missing}}(missing, N)

    # M_prev: total input mass that has already exited before step k₀.
    # cum_m[k₀-1] ≤ M_pipe by definition of k₀, so M_prev = 0 at k₀.
    M_prev = k₀_idx > 1 ? max(0.0, cum_m[k₀_idx - 1] - M_pipe) : 0.0

    for k in k₀_idx:N
        M_new = max(0.0, cum_m[k] - M_pipe)   # input mass exited through end of step k

        i_start = findfirst(x -> x > M_prev, cum_m)
        i_end   = findfirst(x -> x > M_new,  cum_m)
        (isnothing(i_start) || isnothing(i_end)) && break

        dτ₁[k] = k - i_start
        dτ₂[k] = k - i_end

        # α: fraction of plug i_start that exits at step k
        denom_α = ṁ[i_start] * Δt
        α = denom_α > 0 ? (cum_m[i_start] - M_prev) / denom_α : 1.0

        # β: fraction of plug i_end that exits at step k
        cum_m_before_i_end = i_end > 1 ? cum_m[i_end - 1] : 0.0
        denom_β = ṁ[i_end] * Δt
        β = denom_β > 0 ? (M_new - cum_m_before_i_end) / denom_β : 0.0

        αₖ[k] = clamp(α, 0.0, 1.0)
        βₖ[k] = clamp(β, 0.0, 1.0)

        M_prev = M_new
    end

    return (; dτ₁, dτ₂, αₖ, βₖ, k₀=k₀_idx)
end


# For each InsulatedPipe on `path` (producer→load order), compute and return
# the per-pipe delay named tuple from _pipe_time_delays.
function _path_pipe_delays(
    network::Network,
    sim_result::SimulationResults,
    path::Vector{String},
    Δt::Float64,
    N::Int,
)
    delays         = NamedTuple[]
    sim_load_labels = sim_result[:load_labels]

    for i in 1:length(path) - 1
        src, dst = path[i], path[i + 1]
        pipe = network[src, dst]
        pipe isa InsulatedPipe || continue

        r = inner_diameter(pipe) / 2.0
        M = WATER_DENSITY * π * r^2 * pipe_length(pipe)

        # Mass flow through this pipe = sum of all downstream load flows
        ṁ = zeros(Float64, N)
        for dl in downstream_loads(network, dst)
            dl in sim_load_labels || continue
            ṁ .+= sim_result[dl, :mass_flow_load]
        end

        push!(delays, _pipe_time_delays(Δt, M, ṁ))
    end

    return delays
end
