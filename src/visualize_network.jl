# ---------------------------------------------------------------------
# visualization function for Network using GraphMakie
# ---------------------------------------------------------------------

"""
    DEFAULT_ZERO_PIPE_K_ATTRACTION

Default spring constant (dimensionless) for the attraction force pulling a
`ZeroPipe`-connected `LoadNode` toward its source node during auto-positioning.
Forces are normalised by the mean inter-node distance, so this value is
independent of coordinate units.

See also: [`DEFAULT_ZERO_PIPE_K_REPULSION`](@ref), [`compute_zero_pipe_load_positions`](@ref).
"""
const DEFAULT_ZERO_PIPE_K_ATTRACTION = 10.0

"""
    DEFAULT_ZERO_PIPE_K_REPULSION

Default repulsion strength (dimensionless) pushing auto-positioned
`ZeroPipe`-connected `LoadNode`s away from other already-positioned nodes.
Forces are normalised by the mean inter-node distance, so this value is
independent of coordinate units.

See also: [`DEFAULT_ZERO_PIPE_K_ATTRACTION`](@ref), [`compute_zero_pipe_load_positions`](@ref).
"""
const DEFAULT_ZERO_PIPE_K_REPULSION  = 0.3

node_size(::JunctionNode) = 0   # specific size for junction nodes
node_size(::ProducerNode) = 25  # specific size for producer nodes
node_size(::NodeType) = 18      # default size for other node types
node_sizes(mg::MetaGraph) = [node_size(v) for v in vertices_data(mg)]

node_color(::JunctionNode) = colorant"black"
node_color(::ProducerNode) = colorant"green"
node_color(::LoadNode) = colorant"blue"
node_color(::NodeType) = colorant"gray"  # default color
node_colors(mg::MetaGraph) = [node_color(mg[v]) for v in labels(mg)]

node_label(::EmptyNode) = ""  # no label for empty nodes
node_label(n::T) where {T<:NodeType} = n.common.info # label for node types
node_label(::JunctionNode) = ""  # no label for junction nodes
node_label(mg::MetaGraph, i::Int) = node_label(vertex_idx(mg, i))
node_labels(mg::MetaGraph) = [node_label(v) for v in vertices_data(mg)]

edge_colors(mg::MetaGraph) = [edge_color(mg[src, dst]) for (src, dst) in edge_labels(mg)]

edge_linestyle_val(::ZeroPipe) = (:dot, :dense)
edge_linestyle_val(::EdgeType) = :solid
edge_linestyles(mg::MetaGraph) = [edge_linestyle_val(mg[src, dst]) for (src, dst) in edge_labels(mg)]
const max_velocity_for_color = 2.0 # m/s, velocity at which the color will be the most intense
const min_velocity_for_color = 0.0 # m/s, velocity at which the color will be the least intense
const colormap = cgrad(:lajolla) # color gradient for edge coloring based on velocity
function edge_color(e::T)::RGBA{Float64} where {T<:EdgeType} 
    color = colorant"black" # default color for non-pipe edges
    if e isa InsulatedPipe
        vel = water_velocity(e)
        if !ismissing(vel)
            # scale velocity to [1,255//2] for color mapping
            vel_scaled = clamp((vel - min_velocity_for_color) / (max_velocity_for_color - min_velocity_for_color), 0.0, 1.0)
            vel_scaled = 1-vel_scaled # invert so that higher velocity is more intense color
            vel_scaled = Int(round(vel_scaled * (length(colormap.colors)//2 - 1))) + 1 # scale to [1, length(colormap.colors)//2]
            color = colormap[vel_scaled]
        end
    end
    return color
end

function edge_widths(mg::MetaGraph, max_width::Float64=5.0, min_width::Float64=1.0)
    ds = [mg[src,dst] isa InsulatedPipe ? mg[src,dst].physical_params.inner_diameter : missing for (src,dst) in edge_labels(mg)]
    min_d = minimum(ds)
    max_d = maximum(ds)

    if min_d === max_d
        return fill((max_width + min_width) / 2, length(ds)) # if all diameters are the same, return the average width
    end

    # scale edge widths between min_width and max_width based on diameter
    ds_scaled = [ismissing(d) ? min_width : min_width + (d - min_d) / (max_d - min_d) * (max_width - min_width) for d in ds]
    return ds_scaled
end

"""Human-readable edge label used on hover in `visualize_graph!`."""
function edge_info_hover(e::T) where {T<:EdgeType}
    info = e.info * ", L=$(round(pipe_length(e), digits=1)) m, D=$(round(inner_diameter(e)*100, digits=1)) cm"  # label for edge types
    if !ismissing(e.mass_flow)
        info *= ", v=$(round(water_velocity(e), digits=2)) m/s"
    end
    return info
end
edge_info_hover(::EmptyEdge) = ""
edge_info_hover(e::ZeroPipe) = e.info * " (zero pipe)" *
    (!ismissing(e.mass_flow) ? ", ṁ=$(round(e.mass_flow, digits=2)) kg/s" : "")
edge_info(::EmptyEdge) = ""
"""Human-readable edge label used by default in `visualize_graph!`."""
function edge_info(e::T) where {T<:EdgeType} # non-hover label for edge types
    info = ""
    if !ismissing(e.mass_flow)
        info *= "ṁ=$(round(e.mass_flow, digits=2)) kg/s, "
    end
    if !ismissing(e.m_rel)
        info *= "m_rel=$(round(e.m_rel, digits=2))"
    end
    return info
end
"""Vectorized helper returning labels for all edges in a `MetaGraph`."""
edge_infos(mg::MetaGraph) = [edge_info(mg[src, dst]) for (src, dst) in edge_labels(mg)]


"""
    NodeHighlight(; color=nothing, size=nothing)

Specifies visual overrides for a single node in [`highlight_nodes!`](@ref).

Both fields are optional; `nothing` means "keep the default for this node type".

# Fields
- `color`: A Makie-compatible color (e.g. `colorant"red"`, `RGBAf(1,0,0,1)`) or `nothing`.
- `size`: Node marker size in pixels (positive number) or `nothing`.

# Example
```julia
highlights = Dict(
    "Load_1" => NodeHighlight(color=colorant"orange", size=30),
    "Junction_2" => NodeHighlight(size=10),
)
highlight_nodes!(p, nw, highlights)
```
"""
struct NodeHighlight
    color   # nothing or a Makie-compatible color
    size    # nothing or a number
end
NodeHighlight(; color=nothing, size=nothing) = NodeHighlight(color, size)

"""
    highlight_nodes!(p, nw, highlights::Dict{String, NodeHighlight})

Update node colors and sizes in an existing `graphplot` according to `highlights`.

`p` is the third return value of [`visualize_graph!`](@ref).  Each key in
`highlights` is a node label; the associated [`NodeHighlight`](@ref) overrides
its color and/or size.  Nodes not present in `highlights` are unchanged.

# Example
```julia
_, _, p = visualize_graph!(nw)
highlight_nodes!(p, nw, Dict(
    "Load_1" => NodeHighlight(color=colorant"red", size=30),
))
```
"""
function highlight_nodes!(p, nw::Network, highlights::Dict{String, NodeHighlight})
    mg = nw.mg
    colors = copy(p.node_color[])
    sizes  = copy(p.node_size[])
    for (label, h) in highlights
        i = code_for(mg, label)
        if !isnothing(h.color)
            colors[i] = h.color
        end
        if !isnothing(h.size)
            sizes[i] = h.size
        end
    end
    p.node_color[] = colors
    p.node_size[]  = sizes
    return nothing
end

"""
    reset_highlights!(p, nw)

Reset all node colors and sizes in an existing `graphplot` to their type-based
defaults (as computed by `node_colors` and `node_sizes`).

`p` is the third return value of [`visualize_graph!`](@ref).
"""
function reset_highlights!(p, nw::Network)
    mg = nw.mg
    p.node_color[] = node_colors(mg)
    p.node_size[]  = node_sizes(mg)
    return nothing
end

# Returns the point on segment AB closest to point P.
function _closest_point_on_segment(px, py, ax, ay, bx, by)
    dx, dy = bx - ax, by - ay
    len2   = dx^2 + dy^2
    len2 < 1e-10 && return (ax, ay)  # degenerate segment → return A
    t = clamp(((px - ax) * dx + (py - ay) * dy) / len2, 0.0, 1.0)
    return (ax + t * dx, ay + t * dy)
end

"""
    compute_zero_pipe_load_positions(mg; k_attraction, k_repulsion, max_iter) -> Dict{String, Tuple{Float64,Float64}}

Find display positions for `LoadNode`s that have a missing position and are the destination
of a `ZeroPipe` edge whose source is already positioned.

The algorithm places each such load using an attraction-repulsion approach:
- **Attraction** (spring): pulls the load toward its `ZeroPipe` source.
- **Repulsion** (inverse-square): pushes the load away from every other positioned node.

Forces are non-dimensionalised by the mean inter-node distance so that `k_attraction`
and `k_repulsion` behave consistently regardless of coordinate units.

# Keyword arguments
- `k_attraction`: spring constant toward the ZeroPipe source (dimensionless, default [`DEFAULT_ZERO_PIPE_K_ATTRACTION`](@ref)).
- `k_repulsion`: repulsion strength from other nodes (dimensionless, default [`DEFAULT_ZERO_PIPE_K_REPULSION`](@ref)).
- `max_iter`: maximum number of gradient-descent steps per load (default 500).

Returns a `Dict` mapping load label → computed `(x, y)` position. Only loads that
actually need placement are included; nodes that already have a position are not modified.
"""
function compute_zero_pipe_load_positions(mg::MetaGraph;
        k_attraction::Float64 = DEFAULT_ZERO_PIPE_K_ATTRACTION,
        k_repulsion::Float64  = DEFAULT_ZERO_PIPE_K_REPULSION,
        max_iter::Int = 500)

    # Build a working dict of all currently positioned nodes.
    positioned = Dict{String, Tuple{Float64, Float64}}()
    for label in labels(mg)
        p = position(mg[label])
        if !ismissing(p)
            positioned[label] = p
        end
    end

    # Collect loads to place: dst of a ZeroPipe whose src is already positioned.
    to_place = Tuple{String, String}[]
    for (src_lbl, dst_lbl) in edge_labels(mg)
        e = mg[src_lbl, dst_lbl]
        if e isa ZeroPipe && mg[dst_lbl] isa LoadNode && ismissing(mg[dst_lbl].common.position)
            haskey(positioned, src_lbl) && push!(to_place, (src_lbl, dst_lbl))
        end
    end

    isempty(to_place) && return Dict{String, Tuple{Float64, Float64}}()

    # Compute a typical inter-node distance for non-dimensionalisation.
    pos_vals = collect(values(positioned))
    typical = 1.0
    if length(pos_vals) >= 2
        n = length(pos_vals)
        total = 0.0
        for i in 1:n, j in i+1:n
            total += sqrt((pos_vals[i][1] - pos_vals[j][1])^2 + (pos_vals[i][2] - pos_vals[j][2])^2)
        end
        typical = total / (n * (n - 1) / 2)
    end

    # Initialise every load near its source, offset outward from the centre of mass.
    n_fixed = length(positioned)
    cm_x = n_fixed > 0 ? sum(p[1] for p in values(positioned)) / n_fixed : 0.0
    cm_y = n_fixed > 0 ? sum(p[2] for p in values(positioned)) / n_fixed : 0.0

    load_pos = Dict{String, Tuple{Float64, Float64}}()
    for (src_lbl, load_lbl) in to_place
        src = positioned[src_lbl]
        dir_x = src[1] - cm_x
        dir_y = src[2] - cm_y
        d = sqrt(dir_x^2 + dir_y^2)
        if d < 1e-10
            dir_x, dir_y = 1.0, 0.0
        else
            dir_x /= d;  dir_y /= d
        end
        load_pos[load_lbl] = (src[1] + dir_x * typical * 0.3,
                              src[2] + dir_y * typical * 0.3)
    end

    # All positioned fixed nodes — loads repel each other AND these.
    fixed = positioned  # unchanged throughout

    # Collect edge segments between fixed nodes (both endpoints positioned).
    # These are the visible edges in the plot; loads should not land on them.
    edge_segments = Tuple{Tuple{Float64,Float64}, Tuple{Float64,Float64}}[]
    for (src_lbl, dst_lbl) in edge_labels(mg)
        haskey(fixed, src_lbl) && haskey(fixed, dst_lbl) || continue
        push!(edge_segments, (fixed[src_lbl], fixed[dst_lbl]))
    end

    # Converge all loads together in a single shared loop.
    for iter in 1:max_iter
        # Step size decreases over iterations; turbulence amplitude scales the same way.
        α     = typical / (1.0 + iter * 0.02)
        noise = α * 0.05  # turbulence amplitude: 5 % of current step size

        new_pos = copy(load_pos)

        for (src_lbl, load_lbl) in to_place
            src     = positioned[src_lbl]
            pos_x, pos_y = load_pos[load_lbl]
            fx, fy  = 0.0, 0.0

            # Attraction toward ZeroPipe source (linear spring, normalised).
            fx += k_attraction * (src[1] - pos_x) / typical
            fy += k_attraction * (src[2] - pos_y) / typical

            # Repulsion from every fixed node (normalised, inverse-square).
            for (_, other_pos) in fixed
                dx = (pos_x - other_pos[1]) / typical
                dy = (pos_y - other_pos[2]) / typical
                dist2 = dx^2 + dy^2
                dist2 < 1e-6 && continue
                factor = k_repulsion / dist2
                fx += factor * dx
                fy += factor * dy
            end

            # Repulsion from the other loads being placed (using positions from
            # the *previous* iteration so all loads move simultaneously).
            for (other_lbl, other_pos) in load_pos
                other_lbl == load_lbl && continue
                dx = (pos_x - other_pos[1]) / typical
                dy = (pos_y - other_pos[2]) / typical
                dist2 = dx^2 + dy^2
                dist2 < 1e-6 && continue
                factor = k_repulsion / dist2
                fx += factor * dx
                fy += factor * dy
            end

            # Repulsion from the closest point on each fixed edge segment.
            for (a, b) in edge_segments
                cx, cy = _closest_point_on_segment(pos_x, pos_y, a[1], a[2], b[1], b[2])
                dx = (pos_x - cx) / typical
                dy = (pos_y - cy) / typical
                dist2 = dx^2 + dy^2
                dist2 < 1e-6 && continue
                factor = k_repulsion / dist2
                fx += factor * dx
                fy += factor * dy
            end

            # Damped step + small random turbulence to escape local optima.
            new_pos[load_lbl] = (pos_x + α * fx + noise * (rand() - 0.5),
                                 pos_y + α * fy + noise * (rand() - 0.5))
        end

        load_pos = new_pos
    end

    return load_pos
end

"""Visualize a `Network` using GraphMakie.

Returns `(figure, axis, plot)` from `GraphMakie.graphplot`. If edge mass flows
have been computed (e.g. via `steady_state_hydronynamics!`), the plot also
shows flow-dependent edge styling.

`ZeroPipe` edges are drawn with a dotted line style to distinguish them from
physical `InsulatedPipe` edges.

`LoadNode`s that are connected to the network via a `ZeroPipe` and have no
explicit position are automatically placed near their `ZeroPipe` source using
an attraction-repulsion algorithm. Use `k_attraction` and `k_repulsion` to tune
the placement; see [`compute_zero_pipe_load_positions`](@ref) for details.
The module-level constants [`DEFAULT_ZERO_PIPE_K_ATTRACTION`](@ref) and
[`DEFAULT_ZERO_PIPE_K_REPULSION`](@ref) are used as defaults.
"""
function visualize_graph!(nw::Network;
                          k_attraction::Float64 = DEFAULT_ZERO_PIPE_K_ATTRACTION,
                          k_repulsion::Float64  = DEFAULT_ZERO_PIPE_K_REPULSION)
    mg = nw.mg

    # Compute positions for any ZeroPipe-connected loads that lack one.
    auto_pos = compute_zero_pipe_load_positions(mg; k_attraction, k_repulsion)

    node_positions = positions(mg)
    for (label, pos) in auto_pos
        node_positions[code_for(mg, label)] = pos
    end

    f, ax, p = graphplot(mg, layout = any(ismissing.(node_positions)) ? GraphMakie.Spring() : node_positions,
                        node_size = node_sizes(mg),
                        node_color = node_colors(mg),
                        node_attr = (; markerspace = :pixel, marker = :hexagon),
                        nlabels = node_labels(mg),
                        nlabels_attr = (; markerspace = :pixel),
                        nlabels_fontsize = 12,
                        edge_color = edge_colors(mg),
                        edge_width = edge_widths(mg, 10.0, 2.0),
                        edge_attr = (; linestyle = edge_linestyles(mg)),
                        elabels = edge_infos(mg),
                        elabels_fontsize = 12,
                        elabels_attr = (;markerspace = :pixel),
                        arrow_size = 20
                        )
    ax.xautolimitmargin = (0.15, 0.15)
    ax.yautolimitmargin = (0.15, 0.15)
    autolimits!(ax)
    # Equalise x/y data ranges so initial view has 1:1 pixel scale
    fl = ax.finallimits[]
    cx = fl.origin[1] + fl.widths[1] / 2
    cy = fl.origin[2] + fl.widths[2] / 2
    r  = max(fl.widths[1], fl.widths[2]) / 2
    limits!(ax, cx - r, cx + r, cy - r, cy + r)
    hidespines!(ax)
    hidedecorations!(ax)

    # to be closed in callback function
    edge_label_list = [(label_for(nw, src(e)), label_for(nw, dst(e))) for e in edges(nw)]

    function edge_hover_action(state, idx, event, axis)
        # p.edge_width[][idx] = state ? ew[idx] : ew[idx]*1.5
        # p.edge_width[] = p.edge_width[] # trigger observable
        label = edge_label_list[idx]
        
        p.edge_color[][idx] = state ? colorant"red" : edge_color(mg[label...])
        p.edge_color[] = p.edge_color[] # trigger observable
        
        p.elabels[][idx] = state ? edge_info_hover(mg[label...]) : edge_info(mg[label...])
        p.elabels[] = p.elabels[] # trigger observable
    end
    register_interaction!(ax, :ehover, EdgeHoverHandler(edge_hover_action))

    return f, ax, p
end

