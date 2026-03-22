using JLD2

"""
    save_network(filename::AbstractString, network::Network)

Save a `Network` object to a JLD2 file.

The file format is HDF5-based (`.jld2`). All network data — topology, node and edge
types, pipe parameters, plug queues, and load specifications — are serialized.

!!! note "Functions and warnings"
    JLD2 stores functions by module path + name. The built-in functions
    `polynomial_load`, `hockey_load`, and `general_hockey_load` round-trip without
    issues. Anonymous functions and closures defined at the REPL may fail to
    deserialize in a fresh Julia session; prefer named module-level functions.
    JLD2 may print informational warnings about MetaGraphsNext internal anonymous
    functions during save — these are harmless and the file loads correctly.

# Arguments
- `filename`: path to the output file. The `.jld2` extension is appended automatically
  if not already present.
- `network`: the `Network` to save.

# Example
```julia
save_network("my_network.jld2", nw)
```

See also: [`load_network`](@ref).
"""
function save_network(filename::AbstractString, network::Network)
    if !endswith(filename, ".jld2")
        filename = filename * ".jld2"
    end
    jldsave(filename; network)
end

"""
    load_network(filename::AbstractString) -> Network

Load a `Network` object previously saved with [`save_network`](@ref).

# Arguments
- `filename`: path to the `.jld2` file. The `.jld2` extension is appended automatically
  if not already present.

# Example
```julia
nw = load_network("my_network.jld2")
```

See also: [`save_network`](@ref).
"""
function load_network(filename::AbstractString)::Network
    if !endswith(filename, ".jld2")
        filename = filename * ".jld2"
    end
    return load(filename, "network")
end
