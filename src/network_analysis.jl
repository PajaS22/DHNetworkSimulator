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