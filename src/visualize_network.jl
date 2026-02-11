# ---------------------------------------------------------------------
# visualization function for Network using GraphMakie
# ---------------------------------------------------------------------

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

edge_colors(mg::MetaGraph) = [edge_color(mg[src, dst]) for (src, dst) in edge_labels(mg)]  # default edge color
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

    if min_d == max_d
        return fill((max_width + min_width) / 2, length(ds)) # if all diameters are the same, return the average width
    end

    # scale edge widths between min_width and max_width based on diameter
    ds_scaled = [ismissing(d) ? min_width : min_width + (d - min_d) / (max_d - min_d) * (max_width - min_width) for d in ds]
    return ds_scaled
end

function edge_info_hover(e::T) where {T<:EdgeType}
    info = e.info * ", L=$(round(pipe_length(e), digits=1)) m, D=$(round(inner_diameter(e)*100, digits=1)) cm"  # label for edge types
    if !ismissing(e.mass_flow)
        info *= ", v=$(round(water_velocity(e), digits=2)) m/s"
    end
    return info
end
edge_info_hover(::EmptyEdge) = ""
edge_info(::EmptyEdge) = ""
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
edge_infos(mg::MetaGraph) = [edge_info(mg[src, dst]) for (src, dst) in edge_labels(mg)]


function visualize_graph!(nw::Network)
    mg = nw.mg
    # visualize the MetaGraph using GraphMakie
    # it changes the graph properties to store edge ids for interaction handling

    node_positions = positions(mg)

    f, ax, p = graphplot(mg, layout = any(ismissing.(node_positions)) ? GraphMakie.Spring() : node_positions,
                        node_size = node_sizes(mg),
                        node_color = node_colors(mg),
                        node_attr = (; markerspace = :pixel, marker = :hexagon),
                        nlabels = node_labels(mg),
                        nlabels_attr = (; markerspace = :pixel),
                        nlabels_fontsize = 12,
                        edge_color = edge_colors(mg),
                        edge_width = edge_widths(mg, 10.0, 2.0),
                        elabels = edge_infos(mg),
                        elabels_fontsize = 12,
                        elabels_attr = (;markerspace = :pixel),
                        arrow_size = 20
                        )
    ax.aspect = DataAspect()
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

