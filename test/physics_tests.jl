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

    @testset "steady_state_hydrodynamics!" begin
        nw = branching_network()
        total_flow = 15.0
        steady_state_hydrodynamics!(nw, total_flow)

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
        steady_state_hydrodynamics!(nw, 9.0)  # set mass flows first

        velocities = water_velocities(nw)
        @test haskey(velocities, ("producer", "junction"))
        @test haskey(velocities, ("junction", "load1"))
        @test haskey(velocities, ("junction", "load2"))
        # all velocities are positive (flow is moving)
        @test all(v > 0 for v in values(velocities))
    end

    @testset "steady_state_hydrodynamics! with time-varying m_rel" begin
        # Give load1 a 2-step vector: step 1 → m_rel=1, step 2 → m_rel=3
        # At step 1 the split should be 1:2 (as in the constant test).
        # At step 2 the split should be 3:2.
        nw = branching_network()
        set_load_m_rel!(nw, "load1", [1.0, 3.0])
        set_load_m_rel!(nw, "load2", [2.0, 2.0])

        total = 15.0

        steady_state_hydrodynamics!(nw, total, 1)
        @test nw["junction", "load1"].mass_flow ≈ total * (1.0 / 3.0)  atol=1e-10
        @test nw["junction", "load2"].mass_flow ≈ total * (2.0 / 3.0)  atol=1e-10

        steady_state_hydrodynamics!(nw, total, 2)
        @test nw["junction", "load1"].mass_flow ≈ total * (3.0 / 5.0)  atol=1e-10
        @test nw["junction", "load2"].mass_flow ≈ total * (2.0 / 5.0)  atol=1e-10
    end

    @testset "set_m_rel! on pipe edges" begin
        pipe = InsulatedPipe(100.0)
        @test ismissing(pipe.m_rel)
        set_m_rel!(pipe, 1.5)
        @test pipe.m_rel == 1.5

        zpipe = ZeroPipe()
        @test ismissing(zpipe.m_rel)
        set_m_rel!(zpipe, 0.5)
        @test zpipe.m_rel == 0.5

        # vector overload
        vec = [1.0, 2.0, 3.0]
        set_m_rel!(pipe, vec)
        @test pipe.m_rel isa Vector{Float64}
        @test pipe.m_rel == vec
    end

    @testset "m_rel(pipe, step) two-argument accessor" begin
        pipe = InsulatedPipe(100.0)

        # scalar m_rel: step is ignored
        set_m_rel!(pipe, 2.5)
        @test m_rel(pipe, 1)  == 2.5
        @test m_rel(pipe, 99) == 2.5   # step ignored for scalar

        # vector m_rel: indexed by step
        set_m_rel!(pipe, [1.0, 3.0, 5.0])
        @test m_rel(pipe, 1) == 1.0
        @test m_rel(pipe, 2) == 3.0
        @test m_rel(pipe, 3) == 5.0

        zpipe = ZeroPipe()
        set_m_rel!(zpipe, [10.0, 20.0])
        @test m_rel(zpipe, 1) == 10.0
        @test m_rel(zpipe, 2) == 20.0
    end

    @testset "set_relative_mass_flows!(nw) vectorised precompute" begin
        nw = branching_network()
        set_load_m_rel!(nw, "load1", [1.0, 3.0])
        set_load_m_rel!(nw, "load2", [2.0, 2.0])

        set_relative_mass_flows!(nw)

        # Pipes should now carry Vector{Float64} m_rel of length 2
        @test nw["junction", "load1"].m_rel isa Vector{Float64}
        @test nw["junction", "load2"].m_rel isa Vector{Float64}
        @test nw["producer", "junction"].m_rel isa Vector{Float64}
        @test length(nw["junction", "load1"].m_rel) == 2
        @test length(nw["producer", "junction"].m_rel) == 2

        # Values should be correct (normalised so leaf m_rel sums to 1.0 at each step)
        # Step 1: load1=1, load2=2, total=3 → load1=1/3, load2=2/3, junction=1.0
        @test m_rel(nw["junction", "load1"], 1) ≈ 1.0/3.0  atol=1e-10
        @test m_rel(nw["junction", "load2"], 1) ≈ 2.0/3.0  atol=1e-10
        @test m_rel(nw["producer", "junction"], 1) ≈ 1.0    atol=1e-10

        # Step 2: load1=3, load2=2, total=5 → load1=3/5, load2=2/5, junction=1.0
        @test m_rel(nw["junction", "load1"], 2) ≈ 3.0/5.0  atol=1e-10
        @test m_rel(nw["junction", "load2"], 2) ≈ 2.0/5.0  atol=1e-10
        @test m_rel(nw["producer", "junction"], 2) ≈ 1.0    atol=1e-10
    end

    @testset "set_relative_mass_flows!(nw) scalar precompute" begin
        # Constant loads: vectorised form should write scalar Float64 to pipes
        nw = branching_network()  # load1=1.0, load2=2.0

        set_relative_mass_flows!(nw)

        @test nw["junction", "load1"].m_rel isa Float64
        @test nw["junction", "load2"].m_rel isa Float64
        @test nw["producer", "junction"].m_rel isa Float64
        # normalised: load1=1/(1+2)=1/3, load2=2/3, junction=1.0
        @test nw["junction", "load1"].m_rel ≈ 1.0/3.0  atol=1e-10
        @test nw["junction", "load2"].m_rel ≈ 2.0/3.0  atol=1e-10
        @test nw["producer", "junction"].m_rel ≈ 1.0    atol=1e-10
    end

    @testset "Plug fractional k" begin

        # --- collect_exiting_water_plugs!: single full exit ---
        # Plug mass == M_exit → entire plug exits.
        # k_avg_exit = (step-1) + (M_acc + m/2) / M_exit = (5-1) + (0 + 50) / 100 = 4.5
        let plugs = [Plug(80.0, 100.0, 2.0)]
            result = DHN.collect_exiting_water_plugs!(plugs, 100.0, 1.0, 5)
            @test length(result) == 1
            @test result[1].k ≈ 4.5
            @test isempty(plugs)
        end

        # --- collect_exiting_water_plugs!: plug larger than M_exit (split) ---
        # Only 100 of 200 kg exits; remainder keeps its original k.
        let plugs = [Plug(80.0, 200.0, 2.0)]
            result = DHN.collect_exiting_water_plugs!(plugs, 100.0, 1.0, 5)
            @test length(result) == 1
            @test result[1].k ≈ 4.5          # exiting portion: same k_avg as full-exit case
            @test result[1].m ≈ 100.0
            @test length(plugs) == 1
            @test plugs[1].k ≈ 2.0           # remainder preserves original k
            @test plugs[1].m ≈ 100.0
        end

        # --- collect_exiting_water_plugs!: two plugs exit (30 full + 20 partial) ---
        # M_exit=50, plugs: [30kg@k=2, 200kg@k=3]
        # Plug1 (30 kg) exits fully: k_avg = 4 + (0 + 15)/50 = 4.3
        # Plug2 split: 20 kg exits, k_avg = 4 + (30 + 10)/50 = 4.8; 180 kg remains at k=3
        let plugs = [Plug(80.0, 30.0, 2.0), Plug(75.0, 200.0, 3.0)]
            result = DHN.collect_exiting_water_plugs!(plugs, 50.0, 1.0, 5)
            @test length(result) == 2
            @test result[1].k ≈ 4.3
            @test result[2].k ≈ 4.8
            @test length(plugs) == 1
            @test plugs[1].k ≈ 3.0           # remainder of plug2 keeps k=3
            @test plugs[1].m ≈ 180.0
        end

        # --- apply_exit_heat_loss!: continuous transit time ---
        # k_entry=1.0, k_avg_exit=4.5, Δt=1.0 → τ = 3.5 s
        # d=0.1, R=1.0, T_a=10.0
        let d = 0.1, R = 1.0, T_a = 10.0, Δt = 1.0
            τ   = (4.5 - 1.0) * Δt
            A   = π / 4 * d^2
            τ_c = 1000.0 * 4186.0 * A * R
            T_expected = T_a + (80.0 - T_a) * exp(-τ / τ_c)
            plugs = [Plug(80.0, 100.0, 1.0)]
            result = DHN.collect_exiting_water_plugs!(plugs, 100.0, Δt, 5; d=d, R=R, T_a=T_a)
            @test result[1].T ≈ T_expected atol=1e-10
        end

        # --- merge_water_plug_vectors!: k_avg_merged from interval midpoints ---
        # step=3, branch A=[Plug(80,30), Plug(70,30)] (total 60 kg),
        #         branch B=[Plug(40,40)]               (total 40 kg)
        # Change points on f∈[0,1]: {0.5, 1.0}
        # Interval [0,0.5]: A contributes 30 kg @ 80°C, B contributes 20 kg @ 40°C
        #   T = (80*30 + 40*20)/50 = 64.0,  k = (3-1) + (0+0.5)/2 = 2.25
        # Interval [0.5,1.0]: A contributes 30 kg @ 70°C, B contributes 20 kg @ 40°C
        #   T = (70*30 + 40*20)/50 = 58.0,  k = (3-1) + (0.5+1.0)/2 = 2.75
        let pv_A = [Plug(80.0, 30.0, 0.0), Plug(70.0, 30.0, 0.0)],
            pv_B = [Plug(40.0, 40.0, 0.0)]
            result = DHN.merge_water_plug_vectors!([pv_A, pv_B], 3)
            @test length(result) == 2
            @test result[1].T ≈ 64.0
            @test result[1].k ≈ 2.25
            @test result[2].T ≈ 58.0
            @test result[2].k ≈ 2.75
        end

    end  # @testset "Plug fractional k"

    @testset "set_absolute_mass_flows! with step uses m_rel(pipe, step)" begin
        nw = branching_network()
        set_load_m_rel!(nw, "load1", [1.0, 3.0])
        set_load_m_rel!(nw, "load2", [2.0, 2.0])
        set_relative_mass_flows!(nw)  # precompute vectors on pipes

        total = 10.0

        set_absolute_mass_flows!(nw, total, 1)
        @test nw["junction", "load1"].mass_flow ≈ total * (1.0 / 3.0)  atol=1e-10
        @test nw["junction", "load2"].mass_flow ≈ total * (2.0 / 3.0)  atol=1e-10

        set_absolute_mass_flows!(nw, total, 2)
        @test nw["junction", "load1"].mass_flow ≈ total * (3.0 / 5.0)  atol=1e-10
        @test nw["junction", "load2"].mass_flow ≈ total * (2.0 / 5.0)  atol=1e-10
    end

end
