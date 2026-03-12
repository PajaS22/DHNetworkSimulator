using Pkg                                   
Pkg.activate("./scripts")     # activate environment for executing the scripts   
using Revise
using DHNetworkSimulator      # our package in dependency   

using MetaGraphsNext
using Graphs
using GraphMakie, GLMakie
using FileIO
using Plots
using Dates
using TSFrames

Pkg.test("DHNetworkSimulator")            # run tests on DHNetworkSimulator package to make sure everything is working correctly

const DHN = DHNetworkSimulator

network = Network()
junction1_pos = (100.0, 0.0)
junction2_pos = (250.0, 0.0)
network["producer"] = ProducerNode((0.0, 0.0))
network["junction1"] = JunctionNode("J1", junction1_pos)
network["junction2"] = JunctionNode("J2", junction2_pos)
network["load1"] = LoadNode("L1",junction1_pos .+ (25*sqrt(2), 25*sqrt(2)), 2.0)
network["load2"] = LoadNode("L2",junction1_pos .+ (50*sqrt(2), -50*sqrt(2)), 1.0)
network["load3"] = LoadNode("L3",junction2_pos .+ (25*sqrt(2), 25*sqrt(2)), 1.1)
network["load4"] = LoadNode("L4",junction2_pos .+ (50*sqrt(2), -50*sqrt(2)), 3.0)
network["producer", "junction1"] = InsulatedPipe(100)
network["junction1", "junction2"] = InsulatedPipe(150)
network["junction1", "load1"] = InsulatedPipe(50)
network["junction1", "load2"] = InsulatedPipe(100)
network["junction2", "load3"] = InsulatedPipe(50)
network["junction2", "load4"] = InsulatedPipe(100)

network # display the network structure in console

f, ax, p = DHN.visualize_graph!(network)
display(f)

# fill both forward and backward plugs with initial temperatures
fill_pipes_with_initial_temperature!(network, 90.0, 70.0)

# compute steady state flow in all edges for given mass flow from producer
steady_state_hydrodynamics!(network, 10.0)

# ---------------------------------------------------------------------
# TEST SIMULATION 1 - sinusoidal temperature
# ----------------------------------------------------------------------


t = float.(collect(range(0, stop=3*60*60, step=60))) # simulate for three hours with 1 min time step
# sinusoidal mass flow and temperature
function policy(t, Tₐ, T_back)
    mass_flow = 15.0
    temp = 90 + 10*sin(2π*t/(100*60)) # period of 100 minutes, oscillation between 80 and 100 °C
    return ProducerOutput(mass_flow=mass_flow, temperature=temp)
end

results = run_simulation(network, t, policy; T0_f=75.0, T0_b=60.0)

plot_simulation_results(results, :T_load_in; title="Load Inlet Temperatures")    # plot only temperatures
plot_simulation_results(results, ["load1", "load2", "load3", "load4"], :T_load_out; title="Load Outlet Temperatures")  # plot only mass flows
producer_plot = plot_simulation_results(results, :T_producer_in; label="in", title="Producer Temperatures")  # plot mass flows for all loads
plot_simulation_results(producer_plot, results, :T_producer_out; label="out", title="Producer Temperatures")

# ---------------------------------------------------------------------
# TEST SIMULATION 2 - sinusoidal temperature and mass flow
# ----------------------------------------------------------------------

t = float.(collect(range(0, stop=3*60*60, step=60))) # simulate for three hours with 1 min time step
# sinusoidal mass flow and temperature
function policy(t, Tₐ, T_back)
    mass_flow = 15.0 + 5.0*sin(2π*t/(100*60)) # period of 100 minutes, oscillation between 10 and 20 kg/s
    temp = 90 + 10*sin(2π*t/(100*60)) # period of 100 minutes, oscillation between 80 and 100 °C
    return ProducerOutput(mass_flow=mass_flow, temperature=temp)
end


results = run_simulation(network, t, policy; T0_f=75.0, T0_b=60.0)

plot_simulation_results(results, :T_load_in)    # plot only temperatures
plot_simulation_results(results, ["load1", "load2", "load3", "load4"], :T_load_out)  # plot only mass flows
producer_plot = plot_simulation_results(results, :T_producer_in; label="in")  # plot mass flows for all loads
plot_simulation_results(producer_plot, results, :T_producer_out; label="out")
plot_simulation_results(results, :mass_flow_load)

# ---------------------------------------------------------------------
# EXAMPLE 3 - two different load models on individual nodes
# ----------------------------------------------------------------------

# create an artificial outdoor temperature measurement with sinusoidal variation to test the load models
t = float.(collect(range(0, stop=3*60*60, step=60))) # simulate for three hours with 1 min time step
Tₐ = @. -5 + 10*sin(2π*t/(100*60))

function policy(t, Tₐ, T_back)
    return ProducerOutput(mass_flow=40.0, temperature=90.0)
end

# Define a custom exponential load model: high demand at low temperatures, decays exponentially
function exponential_load(params::Vector{Float64}, T_a::Float64)
    # params = [P_max, T_ref, scale]
    # returns kW
    return max(0.0, params[1] * exp(-(T_a - params[2]) / params[3]))
end

# load1 and load2 use the default polynomial model with tweaked coefficients
set_load_fn!(network, "load1", polynomial_load, [600.0, -40.0, 0.8])
set_load_fn!(network, "load2", polynomial_load, [400.0, -20.0, 0.3])
# load3 and load4 use the custom exponential model
set_load_fn!(network, "load3", exponential_load, [500.0, -15.0, 20.0])
set_load_fn!(network, "load4", exponential_load, [800.0, -10.0, 25.0])

results3 = run_simulation(network, t, policy; T0_f=85.0, T0_b=60.0, ambient_temperature=Tₐ)

plot_simulation_results(results3, :power_load; title="Load Power (mixed models)")
plot_simulation_results(results3, :T_load_out; title="Load Outlet Temperatures (mixed models)")

# ---------------------------------------------------------------------
# EXAMPLE 4 - batch set same model with per-load parameters
# ----------------------------------------------------------------------

# Reset all loads to the polynomial model, each with its own calibrated parameters
params_dict = Dict(
    "load1" => [700.0, -45.0, 1.0],
    "load2" => [400.0, -22.0, 0.4],
    "load3" => [420.0, -30.0, 0.6],
    "load4" => [900.0, -50.0, 1.2],
)
set_load_fn!(network, polynomial_load, params_dict)

results4 = run_simulation(network, t, policy; T0_f=90.0, T0_b=60.0, ambient_temperature=Tₐ)

plot_simulation_results(results4, :power_load; title="Load Power (per-load polynomial)")
plot_simulation_results(results4, :T_load_out; title="Load Outlet Temperatures (per-load polynomial)")
plot_simulation_results(results4, :T_producer_in; title="Producer Inlet Temperature (per-load polynomial)")
