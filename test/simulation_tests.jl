@testset "simulation" begin
    network = dhNetwork()
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

    fill_pipes_with_initial_temperature!(network, 90.0, 70.0)
    # check that all pipes are filled with initial temperature
    for e in edges(network)
        @test length(network[src(e), dst(e)].plugs_f) == 1
        @test network[src(e), dst(e)].plugs_f[1].T == 90.0
        @test length(network[src(e), dst(e)].plugs_b) == 1
        @test network[src(e), dst(e)].plugs_b[1].T == 70.0
    end

    t = float.(collect(range(0, stop=3*60*60, step=60))) # simulate for three hours with 1 min time step
    # sinusoidal mass flow and temperature
    function policy(t, Tₐ, T_back)
        mass_flow = 15.0
        temp = 90 + 10*sin(2π*t/(100*60)) # period of 100 minutes, oscillation between 80 and 100 °C
        return ProducerOutput(mass_flow=mass_flow, temperature=temp)
    end
    results = run_simulation(network, t, policy; T0_f=75.0, T0_b=60.0)
    @test length(results) == length(t)
end