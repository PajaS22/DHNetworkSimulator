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

node_size(::JunctionNode) = 0   # junction nodes are invisible
node_size(::ProducerNode) = 25  # pixel size for producer nodes
node_size(::LoadNode) = 18      # pixel size for load nodes
node_size(::SumpNode) = 22      # pixel size for sump nodes
node_size(::NodeType) = 12      # pixel size for other node types
node_sizes(mg::MetaGraph) = [node_size(v) for v in vertices_data(mg)]

node_color(::JunctionNode) = colorant"black"
node_color(::ProducerNode) = colorant"green"
node_color(::LoadNode) = colorant"blue"
node_color(::SumpNode) = colorant"green"
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
        if e.m_rel isa Float64
            info *= "m_rel=$(round(e.m_rel, digits=2))"
        else
            info *= "m_rel=[precomputed, $(length(e.m_rel)) steps]"
        end
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

"""
    EdgeHighlight(; color=nothing, width=nothing)

Specifies visual overrides for a single edge in [`highlight_edges!`](@ref).

Both fields are optional; `nothing` means "keep the default for this edge type".

# Fields
- `color`: A Makie-compatible color (e.g. `colorant"red"`, `RGBAf(1,0,0,1)`) or `nothing`.
- `width`: Edge line width (positive number) or `nothing`.

# Example
```julia
highlights = Dict(
    ("Node_1", "Node_2") => EdgeHighlight(color=colorant"orange", width=5.0),
    ("Node_3", "Node_4") => EdgeHighlight(color=colorant"red"),
)
highlight_edges!(p, nw, highlights)
```
"""
struct EdgeHighlight
    color   # nothing or a Makie-compatible color
    width   # nothing or a number
end
EdgeHighlight(; color=nothing, width=nothing) = EdgeHighlight(color, width)

"""
    highlight_edges!(p, nw, highlights::Dict{Tuple{String,String}, EdgeHighlight})

Update edge colors and widths in an existing `graphplot` according to `highlights`.

`p` is the third return value of [`visualize_graph!`](@ref). Each key in
`highlights` is a `(src_label, dst_label)` tuple; the associated [`EdgeHighlight`](@ref)
overrides its color and/or width. Edges not present in `highlights` are unchanged.

# Example
```julia
_, _, p = visualize_graph!(nw)
highlight_edges!(p, nw, Dict(
    ("Source", "Load_1") => EdgeHighlight(color=colorant"red", width=6.0),
))
```
"""
function highlight_edges!(p, nw::Network, highlights::Dict{Tuple{String,String}, EdgeHighlight})
    mg = nw.mg
    colors = copy(p.edge_color[])
    widths = copy(p.edge_width[])
    for (i, (src_lbl, dst_lbl)) in enumerate(edge_labels(mg))
        key = (src_lbl, dst_lbl)
        haskey(highlights, key) || continue
        h = highlights[key]
        isnothing(h.color) || (colors[i] = h.color)
        isnothing(h.width) || (widths[i] = h.width)
    end
    p.edge_color[] = colors
    p.edge_width[]  = widths
    return nothing
end

"""
    reset_edge_highlights!(p, nw)

Reset all edge colors and widths in an existing `graphplot` to their type-based
defaults (as computed by `edge_colors` and `edge_widths`).

`p` is the third return value of [`visualize_graph!`](@ref).
"""
function reset_edge_highlights!(p, nw::Network)
    mg = nw.mg
    p.edge_color[] = edge_colors(mg)
    p.edge_width[]  = edge_widths(mg, 10.0, 2.0)
    return nothing
end

# Sample-based estimate of the typical inter-node distance (O(1) instead of O(n²)).
function _sample_typical_distance(pos_vals::Vector, n_samples::Int = 200)
    n = length(pos_vals)
    n < 2 && return 1.0
    total = 0.0
    count = 0
    for _ in 1:n_samples
        i, j = rand(1:n), rand(1:n)
        i == j && continue
        total += sqrt((pos_vals[i][1] - pos_vals[j][1])^2 + (pos_vals[i][2] - pos_vals[j][2])^2)
        count += 1
    end
    count == 0 ? 1.0 : total / count
end

"""
    compute_zero_pipe_load_positions(mg; k_attraction, k_repulsion, max_iter, knn_k, knn_seg_k) -> Dict{String, Tuple{Float64,Float64}}

Find display positions for `LoadNode`s that have a missing position and are the destination
of a `ZeroPipe` edge whose source is already positioned.

The algorithm places each such load using an attraction-repulsion approach:
- **Attraction** (spring): pulls the load toward its `ZeroPipe` source.
- **Repulsion** (inverse-square): pushes the load away from the `knn_k` nearest positioned nodes
  and the `knn_seg_k` nearest fixed edge segments.

Forces are non-dimensionalised by the mean inter-node distance so that `k_attraction`
and `k_repulsion` behave consistently regardless of coordinate units.

# Keyword arguments
- `k_attraction`: spring constant toward the ZeroPipe source (dimensionless, default [`DEFAULT_ZERO_PIPE_K_ATTRACTION`](@ref)).
- `k_repulsion`: repulsion strength from other nodes (dimensionless, default [`DEFAULT_ZERO_PIPE_K_REPULSION`](@ref)).
- `max_iter`: maximum number of gradient-descent steps (default 200); stops early when max displacement falls below `1e-3 * typical`.
- `knn_k`: number of nearest fixed nodes used for repulsion per load (default 20; clamped to the number of fixed nodes).
- `knn_seg_k`: number of nearest fixed edge segments used for repulsion per load (default 20; clamped to the number of segments).

Returns a `Dict` mapping load label → computed `(x, y)` position. Only loads that
actually need placement are included; nodes that already have a position are not modified.
"""
function compute_zero_pipe_load_positions(mg::MetaGraph;
        k_attraction::Float64 = DEFAULT_ZERO_PIPE_K_ATTRACTION,
        k_repulsion::Float64  = DEFAULT_ZERO_PIPE_K_REPULSION,
        max_iter::Int  = 200,
        knn_k::Int     = 20,
        knn_seg_k::Int = 20)

    # Build a working dict of all currently positioned nodes.
    positioned = Dict{String, Tuple{Float64, Float64}}()
    for label in labels(mg)
        p = position(mg[label])
        ismissing(p) || (positioned[label] = p)
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

    # Typical distance via random sampling (O(1) instead of O(n²)).
    pos_vals = collect(values(positioned))
    typical  = _sample_typical_distance(pos_vals)

    # Pack fixed positions into a matrix (n_fixed × 2) for fast indexed access.
    fixed_labels_vec   = collect(keys(positioned))
    n_fixed            = length(fixed_labels_vec)
    P                  = Matrix{Float64}(undef, n_fixed, 2)
    for (i, lbl) in enumerate(fixed_labels_vec)
        P[i, 1], P[i, 2] = positioned[lbl]
    end
    fixed_label_to_idx = Dict(lbl => i for (i, lbl) in enumerate(fixed_labels_vec))

    # Centre of mass for initial load placement.
    cm_x = sum(@view P[:, 1]) / n_fixed
    cm_y = sum(@view P[:, 2]) / n_fixed

    # Initialise load positions near their ZeroPipe source, offset outward.
    n_loads     = length(to_place)
    load_labels = [dst for (_, dst) in to_place]
    src_idx_vec = [fixed_label_to_idx[src] for (src, _) in to_place]

    Q     = Matrix{Float64}(undef, n_loads, 2)  # current positions
    Q_new = Matrix{Float64}(undef, n_loads, 2)  # next positions (pre-allocated, no alloc per iter)
    for (i, (src_lbl, _)) in enumerate(to_place)
        sx, sy = positioned[src_lbl]
        dir_x  = sx - cm_x
        dir_y  = sy - cm_y
        d = sqrt(dir_x^2 + dir_y^2)
        if d < 1e-10
            dir_x, dir_y = 1.0, 0.0
        else
            dir_x /= d; dir_y /= d
        end
        Q[i, 1] = sx + dir_x * typical * 0.3
        Q[i, 2] = sy + dir_y * typical * 0.3
    end

    # KNN: for each load find the K nearest fixed nodes once at initialisation.
    # partialsortperm is O(n) and avoids a full sort.
    K       = min(knn_k, n_fixed)
    knn_idx = Vector{Vector{Int}}(undef, n_loads)
    d2_buf  = Vector{Float64}(undef, n_fixed)   # scratch buffer reused below
    for i in 1:n_loads
        @inbounds for j in 1:n_fixed
            d2_buf[j] = (Q[i,1] - P[j,1])^2 + (Q[i,2] - P[j,2])^2
        end
        knn_idx[i] = partialsortperm(d2_buf, 1:K)
    end

    # Precompute edge-segment geometry (fixed throughout):
    # D_seg[j,:] = B[j,:] - A[j,:],  len2_seg[j] = |D_seg[j,:]|²
    segs_src = Tuple{Float64,Float64}[]
    segs_dst = Tuple{Float64,Float64}[]
    for (src_lbl, dst_lbl) in edge_labels(mg)
        haskey(positioned, src_lbl) && haskey(positioned, dst_lbl) || continue
        push!(segs_src, positioned[src_lbl])
        push!(segs_dst, positioned[dst_lbl])
    end
    n_segs   = length(segs_src)
    SegA     = Matrix{Float64}(undef, n_segs, 2)
    D_seg    = Matrix{Float64}(undef, n_segs, 2)
    len2_seg = Vector{Float64}(undef, n_segs)
    for j in 1:n_segs
        SegA[j, 1]  = segs_src[j][1];  SegA[j, 2]  = segs_src[j][2]
        D_seg[j, 1] = segs_dst[j][1] - segs_src[j][1]
        D_seg[j, 2] = segs_dst[j][2] - segs_src[j][2]
        len2_seg[j] = D_seg[j,1]^2 + D_seg[j,2]^2
    end

    # KNN for segments: precompute the K_seg nearest segments (by midpoint) for each load.
    # Avoids iterating all n_segs segments per load per iteration → O(K_seg) instead of O(n_segs).
    K_seg       = min(knn_seg_k, n_segs)
    knn_seg_idx = Vector{Vector{Int}}(undef, n_loads)
    if n_segs > 0
        seg_mid_x  = [SegA[j,1] + D_seg[j,1] * 0.5 for j in 1:n_segs]
        seg_mid_y  = [SegA[j,2] + D_seg[j,2] * 0.5 for j in 1:n_segs]
        seg_d2_buf = Vector{Float64}(undef, n_segs)
        for i in 1:n_loads
            @inbounds for j in 1:n_segs
                seg_d2_buf[j] = (Q[i,1] - seg_mid_x[j])^2 + (Q[i,2] - seg_mid_y[j])^2
            end
            knn_seg_idx[i] = K_seg > 0 ? partialsortperm(seg_d2_buf, 1:K_seg) : Int[]
        end
    else
        for i in 1:n_loads
            knn_seg_idx[i] = Int[]
        end
    end

    # Main optimisation loop.
    # Each load's force update is independent → parallel over loads with Threads.@threads.
    # All threads write to distinct rows of Q_new; reads from Q (previous iteration) are safe.
    for iter in 1:max_iter
        α     = typical / (1.0 + iter * 0.02)
        noise = α * 0.05

        Threads.@threads for i in 1:n_loads
            @inbounds begin
                qx, qy = Q[i, 1], Q[i, 2]
                fx, fy = 0.0, 0.0

                # Attraction toward ZeroPipe source (linear spring, normalised).
                sx = P[src_idx_vec[i], 1];  sy = P[src_idx_vec[i], 2]
                fx += k_attraction * (sx - qx) / typical
                fy += k_attraction * (sy - qy) / typical

                # Repulsion from the K nearest fixed nodes.
                for j in knn_idx[i]
                    dx = (qx - P[j, 1]) / typical
                    dy = (qy - P[j, 2]) / typical
                    d2 = dx*dx + dy*dy
                    d2 < 1e-6 && continue
                    f   = k_repulsion / d2
                    fx += f * dx;  fy += f * dy
                end

                # Repulsion from other loads (previous-iteration positions → simultaneous update).
                for k in 1:n_loads
                    k == i && continue
                    dx = (qx - Q[k, 1]) / typical
                    dy = (qy - Q[k, 2]) / typical
                    d2 = dx*dx + dy*dy
                    d2 < 1e-6 && continue
                    f   = k_repulsion / d2
                    fx += f * dx;  fy += f * dy
                end

                # Repulsion from the closest point on the K_seg nearest fixed edge segments.
                for j in knn_seg_idx[i]
                    len2_seg[j] < 1e-10 && continue
                    t  = clamp(((qx - SegA[j,1]) * D_seg[j,1] + (qy - SegA[j,2]) * D_seg[j,2]) / len2_seg[j], 0.0, 1.0)
                    cx = SegA[j,1] + t * D_seg[j,1]
                    cy = SegA[j,2] + t * D_seg[j,2]
                    dx = (qx - cx) / typical
                    dy = (qy - cy) / typical
                    d2 = dx*dx + dy*dy
                    d2 < 1e-6 && continue
                    f   = k_repulsion / d2
                    fx += f * dx;  fy += f * dy
                end

                Q_new[i, 1] = qx + α * fx + noise * (rand() - 0.5)
                Q_new[i, 2] = qy + α * fy + noise * (rand() - 0.5)
            end
        end

        # Early stopping: terminate when no load moves more than ε (sequential reduction).
        max_disp = 0.0
        @inbounds for i in 1:n_loads
            d = sqrt((Q_new[i,1] - Q[i,1])^2 + (Q_new[i,2] - Q[i,2])^2)
            d > max_disp && (max_disp = d)
        end

        Q, Q_new = Q_new, Q   # swap buffers — zero allocations per iteration

        max_disp < 1e-3 * typical && break
    end

    result = Dict{String, Tuple{Float64, Float64}}()
    for (i, lbl) in enumerate(load_labels)
        result[lbl] = (Q[i, 1], Q[i, 2])
    end
    return result
end

"""
    visualize_graph!(nw; k_attraction, k_repulsion, zoom_factor, zoom_factor_labels) -> (figure, axis, plot)

Visualize a `Network` using GraphMakie.

Returns `(figure, axis, plot)` from `GraphMakie.graphplot`.

**Styling**
- Node markers are fixed-pixel hexagons (unaffected by zoom level).
- `JunctionNode`s remain invisible at all zoom levels.
- `InsulatedPipe` edges are coloured by water velocity when mass flows are
  available; `ZeroPipe` edges are drawn with a dotted line.
- Edge widths are scaled by pipe inner diameter.
- The axis uses `DataAspect()` so that node coordinates are rendered with a
  1:1 aspect ratio throughout zoom and pan. Positions are centred around
  (0, 0) internally to avoid a Makie Float32 precision bug that occurs at
  extreme zoom on large absolute coordinates (e.g. UTM / national grid).

**Zoom-gated labels and arrows**
Node labels, edge labels, and flow arrows are hidden when the visible area is
wide and revealed automatically once the user has zoomed in enough. The
threshold is `zoom_factor × typical_distance`, where `typical_distance` is the
mean inter-node distance sampled from the positioned nodes.

- Default edge label (zoomed in): mass flow and relative mass flow.
- Hover edge label (any zoom): pipe name, length, diameter, and velocity.

**ZeroPipe auto-positioning**
`LoadNode`s connected via a `ZeroPipe` with no explicit position are placed
automatically near their source using an attraction-repulsion algorithm.
See [`compute_zero_pipe_load_positions`](@ref) for tuning options.

# Keyword arguments
- `k_attraction`: spring constant for ZeroPipe load placement (default [`DEFAULT_ZERO_PIPE_K_ATTRACTION`](@ref)).
- `k_repulsion`: repulsion strength for ZeroPipe load placement (default [`DEFAULT_ZERO_PIPE_K_REPULSION`](@ref)).
- `zoom_factor`: node labels and flow arrows appear when the view width drops
  below `zoom_factor × typical_distance` (default `5.0`).
- `zoom_factor_labels`: default edge labels (mass flow, relative mass flow)
  appear when the view width drops below `zoom_factor_labels × typical_distance`
  (default `7.0`). Set higher than `zoom_factor` so edge labels appear before
  arrows as you zoom in. Hover labels are always visible regardless of zoom.
"""
function visualize_graph!(nw::Network;
                          k_attraction::Float64        = DEFAULT_ZERO_PIPE_K_ATTRACTION,
                          k_repulsion::Float64         = DEFAULT_ZERO_PIPE_K_REPULSION,
                          zoom_factor::Float64         = 5.0,
                          zoom_factor_labels::Float64  = 1.0)
    mg = nw.mg

    # Compute positions for any ZeroPipe-connected loads that lack one.
    auto_pos = compute_zero_pipe_load_positions(mg; k_attraction, k_repulsion)

    node_positions = positions(mg)
    for (label, pos) in auto_pos
        node_positions[code_for(mg, label)] = pos
    end

    # Typical inter-node distance — used for the zoom threshold.
    # Computed from raw positions; translation-invariant so the offset below doesn't matter.
    pos_vals = filter(!ismissing, node_positions)
    typical  = _sample_typical_distance(pos_vals)

    # Offset all positions to be centred near (0, 0).
    # Makie's internal Float32 coordinate pipeline (f32c) fails at extreme zoom
    # when absolute coordinate values are large (e.g. UTM / national grid coordinates
    # in the range of hundreds of thousands). Centring removes the large offset while
    # preserving all relative distances and proportions exactly.
    if !isempty(pos_vals)
        x_off = sum(p[1] for p in pos_vals) / length(pos_vals)
        y_off = sum(p[2] for p in pos_vals) / length(pos_vals)
        node_positions = [ismissing(p) ? missing : (p[1] - x_off, p[2] - y_off)
                          for p in node_positions]
    end

    # Start with blank labels and zero arrow size; the zoom callback below
    # reveals them once the user has zoomed in past the threshold.
    n_nodes = nv(mg)
    n_edges = ne(mg)

    f, ax, p = graphplot(mg,
                        layout     = any(ismissing.(node_positions)) ? GraphMakie.Spring() : node_positions,
                        node_size  = node_sizes(mg),
                        node_color = node_colors(mg),
                        node_attr  = (; markerspace = :pixel, marker = :hexagon),
                        nlabels    = fill("", n_nodes),
                        nlabels_attr = (; markerspace = :pixel),
                        nlabels_fontsize = 12,
                        edge_color = edge_colors(mg),
                        edge_width = edge_widths(mg, 10.0, 2.0),
                        edge_attr  = (; linestyle = edge_linestyles(mg)),
                        elabels    = fill("", n_edges),
                        elabels_fontsize = 12,
                        elabels_attr = (; markerspace = :pixel),
                        arrow_size = 0)

    # 1:1 aspect ratio throughout zoom/pan — preserves real-world metre proportions.
    ax.aspect = DataAspect()
    ax.xautolimitmargin = (0.15, 0.15)
    ax.yautolimitmargin = (0.15, 0.15)
    hidespines!(ax)
    hidedecorations!(ax)

    # Precompute node labels and hidden-state arrays (these never change).
    # Edge labels are computed fresh each time the user zooms in so they always
    # reflect the current simulation state, even if run_simulation was called
    # after visualize_graph!.
    lod_node_labels  = node_labels(mg)
    lod_hidden_nodes = fill("", n_nodes)
    lod_hidden_edges = fill("", n_edges)

    # Two independent zoom thresholds:
    #   threshold_arrows — node labels + flow arrows (zoom_factor, default 5×typical)
    #   threshold_labels — default edge labels      (zoom_factor_labels, default 7×typical)
    # Higher factor = appears at a wider view (less zoom required).
    # The early-return guards keep each branch O(1) per frame during smooth zoom.
    threshold_arrows = zoom_factor        * typical
    threshold_labels = zoom_factor_labels * typical

    prev_arrows_visible = Ref{Bool}(false)
    prev_labels_visible = Ref{Bool}(false)

    # apply_lod! does the actual observable updates — called once per settled resize.
    function apply_lod!(w)
        new_arrows = w < threshold_arrows
        new_labels = w < threshold_labels

        if new_arrows != prev_arrows_visible[]
            prev_arrows_visible[] = new_arrows
            p.nlabels[]    = new_arrows ? lod_node_labels : lod_hidden_nodes
            p.arrow_size[] = new_arrows ? 20.0            : 0.0
        end

        if new_labels != prev_labels_visible[]
            prev_labels_visible[] = new_labels
            p.elabels[] = new_labels ? edge_infos(mg) : lod_hidden_edges
        end
    end

    # Debounced LOD observer: during rapid resize (e.g. window maximize) the
    # finallimits observable fires dozens of times per second.  Reacting to every
    # event would queue up many expensive graph re-renders.  Instead we arm a
    # short timer (80 ms) on each event and only call apply_lod! once the window
    # size has settled.
    lod_timer = Ref{Union{Nothing,Timer}}(nothing)

    function update_lod!(limits)
        w = limits.widths[1]
        # Cancel any in-flight timer so we only ever run after the last event.
        t = lod_timer[]
        t !== nothing && close(t)
        lod_timer[] = Timer(_ -> apply_lod!(w), 0.08)
    end
    on(update_lod!, ax.finallimits)

    # autolimits! triggers finallimits → update_lod!; flush immediately so the
    # initial state is set synchronously (no 80 ms delay on first open).
    autolimits!(ax)
    t0 = lod_timer[]
    if t0 !== nothing
        close(t0)
        lod_timer[] = nothing
    end
    apply_lod!(ax.finallimits[].widths[1])

    # Hover: highlight hovered edge in red and show detailed label.
    # Hover labels are always shown on hover regardless of zoom level.
    # On hover-off, restore the default label if the edge-label threshold is met,
    # otherwise restore blank (the LOD callback owns the bulk state).
    edge_label_list = [(label_for(nw, src(e)), label_for(nw, dst(e))) for e in edges(nw)]

    function edge_hover_action(state, idx, event, axis)
        label = edge_label_list[idx]

        # Color: in-place mutation + notify (reliable for this pipeline stage).
        p.edge_color[][idx] = state ? colorant"red" : edge_color(mg[label...])
        notify(p.edge_color)

        # Labels: full array replacement required — the Makie compute pipeline
        # ignores in-place mutations on a Computed value; only setindex! on the
        # input Observable triggers a re-render.
        # Hover labels are shown at any zoom level; on hover-off the default label
        # is restored if the edge-label threshold is met, otherwise blank.
        new_labels = copy(p.elabels[])
        if state
            new_labels[idx] = edge_info_hover(mg[label...])
        else
            within_label_threshold = ax.finallimits[].widths[1] < threshold_labels
            new_labels[idx] = within_label_threshold ? edge_info(mg[label...]) : ""
        end
        p.elabels[] = new_labels
    end
    register_interaction!(ax, :ehover, EdgeHoverHandler(edge_hover_action))

    return f, ax, p
end

