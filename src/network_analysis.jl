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