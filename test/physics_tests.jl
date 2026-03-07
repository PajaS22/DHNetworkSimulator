@testset "Physics" begin

    # Build a simple branching network used across several testsets.
    # Topology: producer → junction → load1 (m_rel=1.0)
    #                               → load2 (m_rel=2.0)
    function branching_network()
        nw = Network()
        nw["producer"] = ProducerNode((0.0, 0.0))
        nw["junction"] = JunctionNode((50.0, 0.0))
        nw["load1"]    = LoadNode("L1", (100.0,  50.0), 1.0)
        nw["load2"]    = LoadNode("L2", (100.0, -50.0), 2.0)
        nw["producer", "junction"] = InsulatedPipe(50)
        nw["junction",  "load1"]   = InsulatedPipe(50)
        nw["junction",  "load2"]   = InsulatedPipe(50)
        return nw
    end

    @testset "power_consumption" begin
        # Default load curve: P(Tₐ) = 540 - 36·Tₐ + 0.6·Tₐ²  [kW]
        load = LoadNode("L", (0.0, 0.0); load=LoadSpec(polynomial_load, [540.0, -36.0, 0.6]))

        # at Tₐ = -10°C: P = (540 + 360 + 60) · 1000 = 960 000 W
        @test power_consumption(load, -10.0) ≈ 960_000.0

        # at Tₐ = 0°C: P = 540 000 W
        @test power_consumption(load, 0.0) ≈ 540_000.0

        # at Tₐ = 30°C: polynomial minimum ≈ 0 (clamped to 0 by max(0, ...))
        @test power_consumption(load, 30.0) ≥ 0.0

        # power_consumption with nothing ambient → 0
        @test power_consumption(load, nothing) == 0.0
    end

    @testset "steady_state_hydronynamics!" begin
        nw = branching_network()
        total_flow = 15.0
        steady_state_hydronynamics!(nw, total_flow)

        # producer → junction carries everything
        @test nw["producer", "junction"].mass_flow ≈ total_flow

        # flow splits proportionally to m_rel (1:2 ratio)
        @test nw["junction", "load1"].mass_flow ≈ total_flow * (1.0 / 3.0)  atol=1e-10
        @test nw["junction", "load2"].mass_flow ≈ total_flow * (2.0 / 3.0)  atol=1e-10

        # flows sum back to total at the junction
        @test (nw["junction", "load1"].mass_flow +
               nw["junction", "load2"].mass_flow) ≈ total_flow  atol=1e-10
    end

    @testset "set_load_params! single node" begin
        nw = branching_network()
        set_load_params!(nw, "load1", 0.7)
        @test nw["load1"].m_rel == 0.7
        @test nw["load2"].m_rel == 2.0  # unchanged
    end

    @testset "set_load_params! dict" begin
        nw = branching_network()
        set_load_params!(nw, Dict("load1" => 0.3, "load2" => 0.8))
        @test nw["load1"].m_rel == 0.3
        @test nw["load2"].m_rel == 0.8
    end

    @testset "water_velocity" begin
        # pipe with no mass flow → missing
        pipe = InsulatedPipe(; length=100.0, inner_diameter=0.1)
        @test ismissing(water_velocity(pipe))

        # pipe with known mass flow: v = ṁ / (ρ · A)
        # A = π·(0.1/2)² ≈ 0.007854 m², ρ = 1000 kg/m³
        pipe2 = InsulatedPipe(; length=100.0, inner_diameter=0.1, mass_flow=1.0)
        A = π * (0.1 / 2)^2
        expected_v = 1.0 / (1000.0 * A)
        @test water_velocity(pipe2) ≈ expected_v  atol=1e-10
    end

    @testset "water_velocities" begin
        nw = branching_network()
        steady_state_hydronynamics!(nw, 9.0)  # set mass flows first

        velocities = water_velocities(nw)
        @test haskey(velocities, ("producer", "junction"))
        @test haskey(velocities, ("junction", "load1"))
        @test haskey(velocities, ("junction", "load2"))
        # all velocities are positive (flow is moving)
        @test all(v > 0 for v in values(velocities))
    end

end
