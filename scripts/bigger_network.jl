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

# const AGG = Aggregation
const DHN = DHNetworkSimulator

# Pkg.test("DHNetworkSimulator")            # run tests on DHNetworkSimulator package to make sure everything is working correctly

digraph, edge_params, node_names, node_positions, power_coefs, m_r = load("datasets/M2_network.jld2", "digraph", "edge_params", "node_names", "node_positions", "power_coefs", "m_r")
# in this file we have temperature measurements from 2024-05-01 to 2025-05-01
outdoor_temperature_ts = load("datasets/outdoor_temperature_ts.jld2", "outdoor_temperature_ts")


network = DHN.Network(digraph)
fill_physical_params!(network, edge_params)
name_nodes!(network, node_names)
identify_producer_and_loads!(network)
fill_node_positions!(network, node_positions)
fill_load_specs!(network, power_coefs, m_r)

f, ax, p = DHN.visualize_graph!(network)
display(f)

# parametrize load nodes by a power consumption function P(T_a), where T_a is the current outdoor temperature, power in kW
# it will be polynomial of 2dn order with specific coefficients for each load node
# P(Tₐ) = p₁ + p₂*Tₐ + p₃*Tₐ²

temperatures = -15:0.5:40
power_plot = Plots.plot(xlabel="Outdoor Temperature (°C)", ylabel="Power Consumption (kW)", title="Power Consumption vs Outdoor Temperature")
for loc in network.load_labels
    P = [power_consumption(network[loc], Tₐ) for Tₐ in temperatures] ./1000.0 # convert from W to kW
    Plots.plot!(power_plot, temperatures, P, label=loc)
end
display(current())

# compute steady state flow in all edges
steady_state_hydronynamics!(network, 100.0)

# ---------------------------------------------------------------------
# TEST SIMULATION 1 - constant input
# ----------------------------------------------------------------------

time_interval = [DateTime(2024, 5, 1), DateTime(2024, 5, 7)]
Tₐ = TSFrames.subset(outdoor_temperature_ts, time_interval...) # select only relevant interval
Tₐ_vec = Tₐ[:, :T_a_avg] # convert TSFrame to vector of outdoor temperatures for each time step

t =  index(Tₐ) # get time vector DateTime format

# constant mass flow and temperature
function policy(t, Tₐ, T_back)
    mass_flow = 110.0
    temp = 80.0
    return ProducerOutput(mass_flow=mass_flow, temperature=temp)
end

results = run_simulation(network, t, policy; T0_f=75.0, T0_b=40.0, ambient_temperature=Tₐ_vec)

# example of accessing results:
# results["M2_VS_1", :mass_flow]
# results["M2_VS_1", :temperature]
# results[:time]
# results[:load_labels]
# results[:load_labels_dict]
# results[:T_producer_in]

# ----------------------------
# example of plotting results:
plot_vs1 = plot_simulation_results(results, ["M2_VS_1"], :T_load_in; label="T_in_VS1", title="Load Temperatures comparison")   # plot only load output temperatures
plot_simulation_results(plot_vs1, results, ["M2_VS_1"], :T_load_out; label="T_out_VS1", title="Load Temperatures comparison")   # add load output temperatures to the same plot
plot_simulation_results(plot_vs1, results, :T_producer_in; lw=2, c=:red, label="Producer T_in", title="Load and Producer Temperatures comparison")    # plot only temperatures
plot_simulation_results(plot_vs1, results, :T_producer_out; lw=2, c=:green, label="Producer T_out", title="Load and Producer Temperatures comparison")    # plot only temperatures

# -----------------------------------------------------------------
# we can also plot whatever we want if we access the result vectors
pwr = results["M2_VS_1", :power_load]
T1 = results["M2_VS_1", :T_load_in]
T2 = results["M2_VS_1", :T_load_out]
Plots.plot([pwr, T1-T2], xlabel="Time Step", ylabel="Power Load (kW)", title="Power Load Over Time")


plot = plot_simulation_results(results, :T_load_out)   # add load output temperatures to the same plot
plot_simulation_results(plot, results, :T_producer_in; lw=2, c=:red)    # plot only temperatures


# ---------------------------------------------------------------------
# TEST SIMULATION 2 - sinusoidal temperature
# ----------------------------------------------------------------------

time_interval = [DateTime(2024, 5, 1), DateTime(2024, 5, 7)]
Tₐ = TSFrames.subset(outdoor_temperature_ts, time_interval...) # select only relevant interval
Tₐ_vec = Tₐ[:, :T_a_avg] # convert TSFrame to vector of outdoor temperatures for each time step

t =  index(Tₐ) # get time vector DateTime format

# sinusoidal mass flow and temperature
function policy(t, Tₐ, T_back)
    t_sec = Dates.value.(t)/1000 # convert from milliseconds to seconds
    mass_flow = 110.0
    temp = 90 + 10*sin(2π*t_sec/(24*3600)) # period of 24 hours, oscillation between 80 and 100 °C
    return ProducerOutput(mass_flow=mass_flow, temperature=temp)
end

results = run_simulation(network, t, policy; T0_f=75.0, T0_b=40.0, ambient_temperature=Tₐ_vec)

# ----------------------------
# example of plotting results:
plot_vs1 = plot_simulation_results(results, ["M2_VS_1"], :T_load_in; label="T_in_VS1", title="Load Temperatures comparison")   # plot only load output temperatures
plot_simulation_results(plot_vs1, results, ["M2_VS_1"], :T_load_out; label="T_out_VS1", title="Load Temperatures comparison")   # add load output temperatures to the same plot
plot_simulation_results(plot_vs1, results, :T_producer_in; lw=2, c=:red, label="Producer T_in", title="Load and Producer Temperatures comparison")    # plot only temperatures
plot_simulation_results(plot_vs1, results, :T_producer_out; lw=2, c=:green, label="Producer T_out", title="Load and Producer Temperatures comparison")    # plot only temperatures

# -----------------------------------------------------------------
plot = plot_simulation_results(results, :T_load_out)   # add load output temperatures to the same plot
plot_simulation_results(plot, results, :T_producer_in; lw=2, c=:red)    # plot only temperatures


# ---------------------------------------------------------------------
# TEST SIMULATION 3 - sinusoidal temperature and triangle mass flow
# ----------------------------------------------------------------------

time_interval = [DateTime(2024, 5, 1), DateTime(2024, 5, 10)]
Tₐ = TSFrames.subset(outdoor_temperature_ts, time_interval...) # select only relevant interval
Tₐ_vec = Tₐ[:, :T_a_avg] # convert TSFrame to vector of outdoor temperatures for each time step

t =  index(Tₐ) # get time vector DateTime format

# sinusoidal mass flow and temperature
function policy(t, Tₐ, T_back)
    t_sec = Dates.value.(t)/1000 # convert from milliseconds to seconds
    mass_flow = 65 + 40 * abs(mod(t_sec / (24*3600), 1) - 0.5)
    temp = 90 + 10*sin(2π*t_sec/(24*3600)) # period of 24 hours, oscillation between 80 and 100 °C
    return ProducerOutput(mass_flow=mass_flow, temperature=temp)
end

results = run_simulation(network, t, policy; T0_f=75.0, T0_b=40.0, ambient_temperature=Tₐ_vec)

plot_simulation_results(results, :mass_flow_producer)   # add load output temperatures to the same plot


# ----------------------------
# example of plotting results:
plot_vs1 = plot_simulation_results(results, ["M2_VS_1"], :T_load_in; label="T_in_VS1", title="Load Temperatures comparison")   # plot only load output temperatures
plot_simulation_results(plot_vs1, results, ["M2_VS_1"], :T_load_out; label="T_out_VS1", title="Load Temperatures comparison")   # add load output temperatures to the same plot
plot_simulation_results(plot_vs1, results, :T_producer_in; lw=2, c=:red, label="Producer T_in", title="Load and Producer Temperatures comparison")    # plot only temperatures
plot_simulation_results(plot_vs1, results, :T_producer_out; lw=2, c=:green, label="Producer T_out", title="Load and Producer Temperatures comparison")    # plot only temperatures

# -----------------------------------------------------------------
plot = plot_simulation_results(results, :T_load_out)   # add load output temperatures to the same plot
plot_simulation_results(plot, results, :T_producer_in; lw=2, c=:red)    # plot only temperatures


# ---------------------------------------------------------------------
# TEST SIMULATION 4 - constant output power with constant mass_flow
# ----------------------------------------------------------------------

time_interval = [DateTime(2024, 5, 1), DateTime(2024, 5, 21)]
Tₐ = TSFrames.subset(outdoor_temperature_ts, time_interval...) # select only relevant interval
Tₐ_vec = Tₐ[:, :T_a_avg] # convert TSFrame to vector of outdoor temperatures for each time step

t =  index(Tₐ) # get time vector DateTime format

producer_power = 4_000_000.0 # 4 MW
c_water = 4186.0 # J/(kg*K)
function policy(t, Tₐ, T_back)
    mass_flow = 80.0
    # P = (T2-T1)*mass_flow*WATER_SPECIFIC_HEAT
    temp = producer_power / (mass_flow * c_water) + T_back
    return ProducerOutput(mass_flow=mass_flow, temperature=temp)
end

results = run_simulation(network, t, policy; T0_f=75.0, T0_b=70.0, ambient_temperature=Tₐ_vec)

# -----------------------------------------------------------------
plot_all = plot_simulation_results(results, :T_load_out; title="Return temperatures on load side")   # add load output temperatures to the same plot
plot_simulation_results(plot_all, results, :T_producer_in; lw=2, c=:red, title="Return temperatures on load side")    # plot only temperatures

# ---------------------------------------------------------------------
plot_power = plot_simulation_results(results, :power_load)   # add load output temperatures to the same plot
plot_power_prod = plot_simulation_results(results, :power_producer; lw=2, c=:red)    # plot only temperatures
pwr_total_loads = sum(results[:power_load], dims=2)[1:end-1]./1000 # sum across loads, exclude last time step which is not valid for producer power, convert from kW to MW
Plots.plot!(plot_power_prod, results[:time][1:end-1], pwr_total_loads; label="Total Power Load")

# ---------------------------------------------------------------------
# TEST SIMULATION 5 - simple temperature control with constant mass flow and maximal power output
# ----------------------------------------------------------------------

time_interval = [DateTime(2024, 5, 1), DateTime(2024, 5, 21)]
Tₐ = TSFrames.subset(outdoor_temperature_ts, time_interval...) # select only relevant interval
Tₐ_vec = Tₐ[:, :T_a_avg] # convert TSFrame to vector of outdoor temperatures for each time step

t =  index(Tₐ) # get time vector DateTime format

T_target = 90
Pwr_max =  5_000_000.0 # 5 MW
function policy(t, Tₐ, T_back)
    t_sec = Dates.value.(t)/1000 # convert from milliseconds to seconds
    mass_flow = 80.0
    P_target = (T_target - T_back) * mass_flow * c_water
    P = min(P_target, Pwr_max)
    temp = P / (mass_flow * c_water) + T_back
    return ProducerOutput(mass_flow=mass_flow, temperature=temp)
end

results = run_simulation(network, t, policy; T0_f=75.0, T0_b=70.0, ambient_temperature=Tₐ_vec)

# ----------------------------
# example of plotting results:
plot_simulation_results(results, :T_producer_out; title="Producer output temperature", linestyle=:solid)   # add load output temperatures to the same plot
plot_vs1 = plot_simulation_results(results, ["M2_VS_1"], :T_load_in; label="T_in_VS1", title="Load Temperatures comparison")   # plot only load output temperatures
plot_simulation_results(plot_vs1, results, ["M2_VS_1"], :T_load_out; label="T_out_VS1", title="Load Temperatures comparison")   # add load output temperatures to the same plot
plot_simulation_results(plot_vs1, results, :T_producer_in; lw=2, c=:red, label="Producer T_in", title="Load and Producer Temperatures comparison")    # plot only temperatures
plot_simulation_results(plot_vs1, results, :T_producer_out; lw=2, c=:green, label="Producer T_out", title="Load and Producer Temperatures comparison")    # plot only temperatures

# -----------------------------------------------------------------
plot_all = plot_simulation_results(results, :T_load_out; title="Return temperatures")   # add load output temperatures to the same plot
plot_simulation_results(plot_all, results, :T_producer_in; title="Return temperatures", lw=2, c=:red)    # plot only temperatures

# ---------------------------------------------------------------------
plot_power = plot_simulation_results(results, :power_load)   # add load output temperatures to the same plot
plot_power_prod = plot_simulation_results(results, :power_producer; lw=2, c=:red)    # plot only temperatures
pwr_total_loads = sum(results[:power_load], dims=2)[1:end-1]./1000 # sum across loads, exclude last time step which is not valid for producer power, convert from kW to MW
Plots.plot!(plot_power_prod, results[:time][1:end-1], pwr_total_loads; label="Total Power Load")
Energy_producer = sum(results[:power_producer][1:end-1]) * 600 / 3600.0 # in MWh, sum across time steps, convert from seconds to hours
Energy_loads = sum(pwr_total_loads) * 600 / 3600.0 # in MWh, sum across time steps, convert from seconds to hours
println("Total energy produced: $(round(Energy_producer, digits=2)) MWh")
println("Total energy consumed by loads: $(round(Energy_loads, digits=2)) MWh")
println("Energy balance (Produced - Consumed): $(round(Energy_producer - Energy_loads, digits=2)) MWh")
