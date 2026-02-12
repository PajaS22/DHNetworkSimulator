"""Plot time series from [`SimulationResults`](@ref).

`plot_simulation_results` is a small helper around Plots.jl that plots one physical variable either:

- for selected **load** labels, or
- for the **producer** (for producer-only variables).

```julia
plot_simulation_results(sr::SimulationResults, physical_var::Symbol; kwargs...)
plot_simulation_results(sr::SimulationResults, labels::Vector{String}, physical_var::Symbol; kwargs...)
plot_simulation_results(plot::Plots.Plot{Plots.GRBackend}, sr::SimulationResults, physical_var::Symbol; kwargs...)
plot_simulation_results(plot::Plots.Plot{Plots.GRBackend}, sr::SimulationResults, labels::Vector{String}, physical_var::Symbol; kwargs...)
```

# Arguments
- `sr::SimulationResults`: simulation output from [`run_simulation`](@ref).
- `labels::Vector{String}`: load labels to include in plot (only used for load variables),
                            if not provided, all load labels are plotted.
- `physical_var::Symbol`: what values to plot.
    Supported options are:
    - temperatures: `:T_load_in`, `:T_load_out`, `:T_producer_in`, `:T_producer_out`
    - mass flows: `:mass_flow_load`, `:mass_flow_producer`
    - powers: `:power_load`, `:power_producer`
- `plot::Plots.Plot{Plots.GRBackend}`: an existing plot to add lines to,
                                       if not provided, a new plot is created.

# Keyword Arguments
- `kwargs...`: forwarded to `Plots.plot!` (e.g. `linewidth`, `color`, `legend`, ...).

# Returns
- The `Plots.Plot` object.

# Notes
- If `sr[:time]` is `Vector{Float64}`, the x-axis is converted from seconds to minutes.
- `:power_producer` has length `N-1` and is plotted against `time[1:end-1]`.
- When plotting producer variables, the line style defaults to dashed (unless you override `linestyle`).

# Examples
```julia
sr = run_simulation(network, t, policy)

plot_simulation_results(sr, :T_load_in)
plot_simulation_results(sr, ["L1", "L2"], :mass_flow_load)
plot_simulation_results(sr, :power_producer)
```
"""
function plot_simulation_results(plot::Plots.Plot{Plots.GRBackend}, sr::SimulationResults, labels::Vector{String}, physical_var::Symbol; kwargs...)
    # Plot results for specified load node labels and physical variable
    # options for physical_var: :T_load_in, :T_load_out, :mass_flow, :power
    # if flag "producer" is true, also plot the input variable (T/mass_flow/power) on the side of producer

    options = (:T_load_in, :T_load_out, :T_producer_in, :T_producer_out, :mass_flow_load, :mass_flow_producer, :power_load, :power_producer)
    if !(physical_var ∈ options)
        error("Invalid physical_var. Valid options are: $(options)")
    end

    if(physical_var ∈ (:T_load_in, :T_load_out, :T_producer_in, :T_producer_out))
        Plots.ylabel!(plot,"Temperature (°C)")
        Plots.title!(plot, "Temperatures Over Time")
    end
    if(physical_var ∈ (:mass_flow_load, :mass_flow_producer))
        Plots.ylabel!(plot, "Mass Flow (kg/s)")
        Plots.title!(plot, "Mass Flows Over Time")
    end
    if(physical_var == :power_load)
        Plots.ylabel!(plot, "Power Consumption (kW)")
        Plots.title!(plot, "Power Consumption Over Time")
    end
    if(physical_var == :power_producer)
        Plots.ylabel!(plot, "Power Consumption (MW)")
        Plots.title!(plot, "Power Consumption Over Time")
    end
    
    # check time vector format
    time = sr[:time] # time in seconds
    if(eltype(time) <: DateTime)
        time_plot = time # dont convert
    elseif eltype(time) <: Float64
        time_plot = time ./ 60 # convert to time in minutes, assume time is initially in seconds
        Plots.xlabel!(plot, "Time (min)")
    else
        error("Unsupported time vector element type: $(eltype(time)). Expected Float64 or DateTime.")
    end
    
    producer_vars = (:T_producer_in, :T_producer_out, :mass_flow_producer, :power_producer)
    load_vars = (:T_load_in, :T_load_out, :mass_flow_load, :power_load)

    if (physical_var ∈ producer_vars)
        if :linestyle ∉ keys(kwargs)
            kwargs = merge(NamedTuple(kwargs), (linestyle=:dash,))
        end
        if physical_var == :power_producer
            Plots.plot!(plot, time_plot[1:end-1], sr[:power_producer]; label="producer", kwargs...)
        else
            Plots.plot!(plot, time_plot, sr[physical_var]; label="producer", kwargs...)
        end
    elseif (physical_var ∈ load_vars)
        for label in labels
            Plots.plot!(plot, time_plot, sr[label, physical_var]; label=label, linestyle=:solid, kwargs...)
        end
    else
        error("Invalid physical_var. Valid options are: $(options)")
    end
    return plot
end

# Plot all load nodes
plot_simulation_results(plot::Plots.Plot{Plots.GRBackend}, sr::SimulationResults, physical_var::Symbol; kwargs...) = plot_simulation_results(plot, sr, collect(sr[:load_labels]), physical_var; kwargs...)
plot_simulation_results(sr::SimulationResults, physical_var::Symbol; kwargs...) = plot_simulation_results(sr, collect(sr[:load_labels]), physical_var; kwargs...)

# initialize plot object
function plot_simulation_results(sr::SimulationResults, labels::Vector{String}, physical_var::Symbol; kwargs...)
    # initialize plot object
    plot = Plots.plot(xlabel="Time", ylabel="", title="", legend=:outerright; kwargs...)
    return plot_simulation_results(plot, sr, labels, physical_var; kwargs...)
end

