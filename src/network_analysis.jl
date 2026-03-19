"""
    producer_loads_volumes(network::Network) -> Dict{String, Float64}

Return the total pipe volume (m³) along the path from the producer to each
load node in the network.

Each entry maps a load label to the sum of pipe volumes on the unique path
from the producer to that load. This is useful for estimating delay times or
the amount of water in transit to each consumer.

# Example
```julia
vols = producer_loads_volumes(network)
# => Dict("load_1" => 0.045, "load_2" => 0.12, ...)
```
"""
function producer_loads_volumes(network::Network)::Dict{String, Float64}
    loads = collect(network.load_labels)
    mg = network.mg
    producer_id = index_for(network, network.producer_label)
    D = zeros(Float64, nv(network), nv(network))
    # fill in the volumes of the pipes in the adjacency matrix
    for e in edges(network)
        p = network[src(e), dst(e)]
        D[src(e), dst(e)] = volume(p)
    end
    ds = dijkstra_shortest_paths(mg, producer_id, D);
    load_dists = [ds.dists[index_for(network, load_label)] for load_label in loads]
    return Dict(loads[i] => load_dists[i] for i in 1:length(loads))
end

"""
    producer_loads_delays(network::Network) -> Dict{String, Float64}

Return the true plug-flow transit delay (seconds) from the producer to each load,
accounting for the actual mass flow through every pipe segment.

Unlike [`producer_loads_volumes`](@ref), which sums raw pipe volumes, this function
weights each segment by its own transit time `V/ṁ`.  In a chain or branching network
the mass flow decreases towards the leaves (each tap removes flow), so the true delay
to a distant load is significantly larger than `ΣV / ṁ_total`.

Requires `steady_state_hydrodynamics!` (or equivalent) to have been called so that
every `InsulatedPipe` has a non-missing, positive `mass_flow`.  `ZeroPipe` edges
contribute zero delay (they are instantaneous).

# Example
```julia
steady_state_hydrodynamics!(network, 100.0)
delays = producer_loads_delays(network)
max_delay_min = maximum(values(delays)) / 60
println("Longest transit delay: \$(round(max_delay_min, digits=1)) min")
```
"""
function producer_loads_delays(network::Network)::Dict{String, Float64}
    loads = collect(network.load_labels)
    mg = network.mg
    producer_id = index_for(network, network.producer_label)
    D = zeros(Float64, nv(network), nv(network))
    for e in edges(network)
        p = network[src(e), dst(e)]
        if p isa InsulatedPipe
            mf = mass_flow(p)
            # τ = V·ρ / ṁ  [s].  WATER_DENSITY converts m³ to kg so that τ [s] = V[m³]·ρ[kg/m³] / ṁ[kg/s].
            D[src(e), dst(e)] = (!ismissing(mf) && mf > 0.0) ? volume(p) * WATER_DENSITY / mf : 0.0
        end
        # ZeroPipe → delay 0 (already zero from zeros())
    end
    ds = dijkstra_shortest_paths(mg, producer_id, D)
    load_dists = [ds.dists[index_for(network, l)] for l in loads]
    return Dict(loads[i] => load_dists[i] for i in 1:length(loads))
end