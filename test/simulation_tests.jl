@testset "simulation" begin

    # ─────────────────────────────────────────────────────────────────────────
    # Shared fixture
    # ─────────────────────────────────────────────────────────────────────────

    # Simple 2-load branching network used by most testsets.
    # Topology: producer → (pipe 50 m) → junction → (pipe 50 m) → load1 (m_rel=1.0)
    #                                               → (pipe 50 m) → load2 (m_rel=2.0)
    #
    # The unequal m_rel values (1.0 vs 2.0) ensure the hydraulic split is non-trivial:
    # load1 receives 1/3 of total flow, load2 receives 2/3.  Equal pipes (all 50 m) keep
    # transit times the same across branches so transport-delay tests stay simple.
    function make_network()
        nw = Network()
        nw["producer"] = ProducerNode((0.0, 0.0))           # heat source at origin
        nw["junction"] = JunctionNode((50.0, 0.0))          # single branch point, 50 m from producer
        nw["load1"]    = LoadNode("L1", (100.0,  50.0), 1.0)  # m_rel=1.0 → 1/3 of total flow
        nw["load2"]    = LoadNode("L2", (100.0, -50.0), 2.0)  # m_rel=2.0 → 2/3 of total flow
        nw["producer", "junction"] = InsulatedPipe(50)      # supply main, 50 m
        nw["junction",  "load1"]   = InsulatedPipe(50)      # branch to load1, 50 m
        nw["junction",  "load2"]   = InsulatedPipe(50)      # branch to load2, 50 m
        return nw
    end

    t  = float.(collect(range(0, stop=60*60, step=60)))  # 1 h simulation, 1-min time steps (N=61)
    N  = length(t)
    Tₐ = fill(5.0, N)   # constant 5 °C outdoor temperature — cold enough for nonzero load power

    # Constant return temperatures injected at each load for :backward_only and :hybrid tests.
    # Values are typical DH return temps (35–40 °C) and are lower than the 60–80 °C supply,
    # so hybrid power (ṁ·Cₚ·(T_in − T_out)) is positive under normal conditions.
    T_inject = Dict("load1" => fill(40.0, N), "load2" => fill(35.0, N))

    # ─────────────────────────────────────────────────────────────────────────
    # ProducerOutput struct
    # ─────────────────────────────────────────────────────────────────────────

    @testset "ProducerOutput struct" begin
        # positional constructor
        p = ProducerOutput(15.0, 90.0)
        @test p.mass_flow   == 15.0
        @test p.temperature == 90.0

        # keyword constructor with explicit temperature
        p2 = ProducerOutput(mass_flow=10.0, temperature=80.0)
        @test p2.mass_flow   == 10.0
        @test p2.temperature == 80.0

        # keyword constructor — temperature defaults to nothing
        p3 = ProducerOutput(mass_flow=10.0)
        @test p3.temperature === nothing

        # explicit nothing
        p4 = ProducerOutput(mass_flow=10.0, temperature=nothing)
        @test p4.temperature === nothing

        # nothing is accepted in the Union type
        @test p4.temperature isa Nothing
    end

    # ─────────────────────────────────────────────────────────────────────────
    # Input validation
    # ─────────────────────────────────────────────────────────────────────────

    @testset "run_simulation validation" begin
        nw      = make_network()
        p_full  = (t, Ta, Tb) -> ProducerOutput(mass_flow=10.0, temperature=80.0)
        p_2arg  = (t, Ta, Tb) -> ProducerOutput(mass_flow=10.0, temperature=80.0)

        # invalid mode symbol
        @test_throws ErrorException run_simulation(nw, t, p_full; mode=:invalid)

        # :backward_only and :hybrid require T_return_inject
        @test_throws ErrorException run_simulation(nw, t, p_2arg; mode=:backward_only)
        @test_throws ErrorException run_simulation(nw, t, p_2arg; mode=:hybrid)

        # T_return_inject with a key that is not a load label
        @test_throws ErrorException run_simulation(nw, t, p_2arg;
            mode=:backward_only,
            T_return_inject=Dict("nonexistent" => fill(40.0, N), "load2" => fill(35.0, N)))

        # T_return_inject vector with wrong length
        @test_throws ErrorException run_simulation(nw, t, p_2arg;
            mode=:backward_only,
            T_return_inject=Dict("load1" => fill(40.0, N - 1), "load2" => fill(35.0, N)))

        # policy Vector with wrong length
        short_vec = [ProducerOutput(10.0, 80.0) for _ in 1:N-1]
        @test_throws ErrorException run_simulation(nw, t, short_vec; mode=:full)

        # policy Vector with temperature=nothing in :full
        nil_vec = [ProducerOutput(10.0, nothing) for _ in 1:N]
        @test_throws ErrorException run_simulation(nw, t, nil_vec; mode=:full)

        # policy Vector with temperature=nothing in :forward_only
        @test_throws ErrorException run_simulation(nw, t, nil_vec; mode=:forward_only)

        # policy Vector with temperature=nothing in :hybrid
        @test_throws ErrorException run_simulation(nw, t, nil_vec;
            mode=:hybrid, T_return_inject=T_inject)
    end

    # ─────────────────────────────────────────────────────────────────────────
    # :full mode — structure and basic physics
    # ─────────────────────────────────────────────────────────────────────────

    @testset ":full mode" begin
        nw = make_network()
        p  = (t, Ta, Tb) -> ProducerOutput(mass_flow=10.0, temperature=80.0)
        sr = run_simulation(nw, t, p; mode=:full, T0_f=60.0, T0_b=50.0,
                            ambient_temperature=Tₐ)

        @test length(sr) == N

        # result matrix shapes
        @test size(sr.mass_flow_load)  == (N, 2)
        @test size(sr.T_load_in)       == (N, 2)
        @test size(sr.T_load_out)      == (N, 2)
        @test size(sr.power_load)      == (N, 2)
        @test length(sr.T_producer_in)  == N
        @test length(sr.T_producer_out) == N
        @test length(sr.power_producer) == N - 1

        # no optional field is Nothing in :full mode
        @test !isnothing(sr.T_load_out)
        @test !isnothing(sr.T_producer_in)
        @test !isnothing(sr.power_load)
        @test !isnothing(sr.power_producer)

        # no NaN in populated fields
        @test !any(isnan, sr.T_load_in)
        @test !any(isnan, sr.T_load_out)
        @test !any(isnan, sr.T_producer_in)

        # producer output temperature matches policy exactly
        @test all(sr.T_producer_out .≈ 80.0)

        # supply temperatures are in a physical range
        @test all(55.0 .< sr.T_load_in .< 85.0)

        # return temperatures are below supply temperatures
        @test all(sr.T_load_out .< sr.T_load_in)

        # load power is non-negative (clamped by power_consumption)
        @test all(sr.power_load .>= 0.0)

        # load labels indexing
        @test "load1" ∈ keys(sr[:load_labels_dict])
        @test "load2" ∈ keys(sr[:load_labels_dict])
        @test sr["load1", :T_load_in] isa Vector{Float64}
        @test length(sr["load1", :T_load_in]) == N
    end

    # ─────────────────────────────────────────────────────────────────────────
    # :full mode with Vector{ProducerOutput} policy
    # ─────────────────────────────────────────────────────────────────────────

    @testset ":full mode - Vector policy" begin
        nw   = make_network()
        pvec = [ProducerOutput(10.0, 80.0) for _ in 1:N]
        sr   = run_simulation(nw, t, pvec; mode=:full, ambient_temperature=Tₐ)

        @test length(sr) == N
        @test !isnothing(sr.T_producer_in)
        @test !isnothing(sr.power_load)
        @test !isnothing(sr.power_producer)
        @test all(sr.T_producer_out .≈ 80.0)
    end

    # ─────────────────────────────────────────────────────────────────────────
    # :forward_only mode
    # ─────────────────────────────────────────────────────────────────────────

    @testset ":forward_only mode" begin
        nw = make_network()
        p  = (t, Ta, Tb) -> ProducerOutput(mass_flow=10.0, temperature=80.0)
        sr = run_simulation(nw, t, p; mode=:forward_only, T0_f=60.0,
                            ambient_temperature=Tₐ)

        @test length(sr) == N
        @test size(sr.T_load_in) == (N, 2)

        # supply-side fields are populated and finite
        @test !any(isnan, sr.T_load_in)
        @test all(sr.T_producer_out .≈ 80.0)
        @test all(sr.mass_flow_load .> 0.0)

        # backward fields are absent
        @test isnothing(sr.T_load_out)
        @test isnothing(sr.T_producer_in)
        @test isnothing(sr.power_load)
        @test isnothing(sr.power_producer)

        # :forward_only also works with a Vector policy
        pvec = [ProducerOutput(10.0, 80.0) for _ in 1:N]
        sr2  = run_simulation(nw, t, pvec; mode=:forward_only)
        @test isnothing(sr2.T_producer_in)
        @test !any(isnan, sr2.T_load_in)
    end

    # ─────────────────────────────────────────────────────────────────────────
    # :backward_only mode — temperature=nothing (no supply)
    # ─────────────────────────────────────────────────────────────────────────

    @testset ":backward_only mode - temperature=nothing" begin
        nw   = make_network()
        pvec = [ProducerOutput(10.0, nothing) for _ in 1:N]
        sr   = run_simulation(nw, t, pvec; mode=:backward_only,
                              T_return_inject=T_inject, T0_b=50.0,
                              ambient_temperature=Tₐ)

        @test length(sr) == N

        # T_load_in is all NaN (forward step skipped)
        @test all(isnan, sr.T_load_in)

        # T_producer_out is all NaN (no supply temperature)
        @test all(isnan, sr.T_producer_out)

        # power_load is absent
        @test isnothing(sr.power_load)

        # power_producer is nothing (T_producer_out is all NaN)
        @test isnothing(sr.power_producer)

        # T_load_out equals the injected values exactly
        col1 = sr[:load_labels_dict]["load1"]
        col2 = sr[:load_labels_dict]["load2"]
        @test all(sr.T_load_out[:, col1] .≈ 40.0)
        @test all(sr.T_load_out[:, col2] .≈ 35.0)

        # T_producer_in is populated and finite (backward step ran)
        @test !isnothing(sr.T_producer_in)
        @test !any(isnan, sr.T_producer_in)

        # mass flows are positive (hydraulics ran)
        @test all(sr.mass_flow_load .> 0.0)

        # :backward_only also accepts a function policy (T_back is missing in this mode)
        p_fn = (t, Ta, Tb) -> ProducerOutput(mass_flow=10.0, temperature=nothing)
        sr2  = run_simulation(nw, t, p_fn; mode=:backward_only,
                              T_return_inject=T_inject, T0_b=50.0)
        @test !isnothing(sr2.T_producer_in)
    end

    # ─────────────────────────────────────────────────────────────────────────
    # :backward_only mode — with non-nothing temperature → power_producer computed
    # ─────────────────────────────────────────────────────────────────────────

    @testset ":backward_only mode - temperature provided" begin
        nw   = make_network()
        pvec = [ProducerOutput(10.0, 80.0) for _ in 1:N]
        sr   = run_simulation(nw, t, pvec; mode=:backward_only,
                              T_return_inject=T_inject, T0_b=50.0)

        @test all(sr.T_producer_out .≈ 80.0)
        @test !isnothing(sr.power_producer)
        @test length(sr.power_producer) == N - 1
        # T_load_in is still NaN (forward step skipped regardless)
        @test all(isnan, sr.T_load_in)
    end

    # ─────────────────────────────────────────────────────────────────────────
    # :hybrid mode
    # ─────────────────────────────────────────────────────────────────────────

    @testset ":hybrid mode" begin
        nw = make_network()
        p  = (t, Ta, Tb) -> ProducerOutput(mass_flow=10.0, temperature=80.0)
        sr = run_simulation(nw, t, p; mode=:hybrid, T_return_inject=T_inject,
                            T0_f=60.0, T0_b=50.0, ambient_temperature=Tₐ)

        @test length(sr) == N

        # forward fields populated
        @test !any(isnan, sr.T_load_in)
        @test all(sr.T_producer_out .≈ 80.0)

        # backward fields populated
        @test !isnothing(sr.T_load_out)
        @test !isnothing(sr.T_producer_in)
        @test !any(isnan, sr.T_producer_in)
        @test !isnothing(sr.power_load)
        @test !isnothing(sr.power_producer)

        # T_load_out equals injected values exactly
        col1 = sr[:load_labels_dict]["load1"]
        col2 = sr[:load_labels_dict]["load2"]
        @test all(sr.T_load_out[:, col1] .≈ 40.0)
        @test all(sr.T_load_out[:, col2] .≈ 35.0)

        # power_load formula: P [kW] = ṁ · Cₚ · (T_in − T_out) / 1000
        Cp = DHNetworkSimulator.WATER_SPECIFIC_HEAT
        P_expected = sr.mass_flow_load[:, col1] .* Cp .*
                     (sr.T_load_in[:, col1] .- sr.T_load_out[:, col1]) ./ 1000.0
        @test sr.power_load[:, col1] ≈ P_expected  atol=1e-8

        # also works with a Vector policy
        pvec = [ProducerOutput(10.0, 80.0) for _ in 1:N]
        sr2  = run_simulation(nw, t, pvec; mode=:hybrid, T_return_inject=T_inject,
                              T0_f=60.0, T0_b=50.0)
        @test !isnothing(sr2.power_load)
        @test all(sr2.T_load_out[:, col1] .≈ 40.0)
    end

    # ─────────────────────────────────────────────────────────────────────────
    # :hybrid — negative power is allowed (not clamped)
    # ─────────────────────────────────────────────────────────────────────────

    @testset ":hybrid - negative power preserved" begin
        nw = make_network()
        # Inject return temperature higher than initial forward fill temperature (60 °C).
        # At early steps T_load_in ≈ 60 °C < 90 °C inject → P = m·Cp·(60−90)/1000 < 0.
        T_hot = Dict("load1" => fill(90.0, N), "load2" => fill(90.0, N))
        p  = (t, Ta, Tb) -> ProducerOutput(mass_flow=10.0, temperature=80.0)
        sr = run_simulation(nw, t, p; mode=:hybrid, T_return_inject=T_hot,
                            T0_f=60.0, T0_b=50.0)

        col1 = sr[:load_labels_dict]["load1"]
        # At step 1, T_load_in ≈ T0_f = 60 < 90 → power should be negative
        @test sr.power_load[1, col1] < 0.0
        # In general, negative values exist and are not clamped to zero
        @test any(sr.power_load .< 0.0)
    end

    # ─────────────────────────────────────────────────────────────────────────
    # Chaining: hybrid → backward_only
    # ─────────────────────────────────────────────────────────────────────────

    @testset "chaining hybrid → backward_only" begin
        nw = make_network()
        p  = (t, Ta, Tb) -> ProducerOutput(mass_flow=10.0, temperature=80.0)
        sr_hyb = run_simulation(nw, t, p; mode=:hybrid, T_return_inject=T_inject,
                                T0_f=60.0, T0_b=50.0)

        # Extract T_load_out from hybrid run as inject for a backward-only follow-up
        T_inj2 = Dict(l => sr_hyb[l, :T_load_out]
                      for l in keys(sr_hyb[:load_labels_dict]))

        p_bwd  = (t, Ta, Tb) -> ProducerOutput(mass_flow=10.0)
        sr_bwd = run_simulation(nw, t, p_bwd; mode=:backward_only,
                                T_return_inject=T_inj2, T0_b=50.0)

        @test length(sr_bwd) == N
        @test !isnothing(sr_bwd.T_producer_in)
        @test !any(isnan, sr_bwd.T_producer_in)

        # T_load_out in backward run equals T_load_out from hybrid
        col1 = sr_bwd[:load_labels_dict]["load1"]
        @test sr_bwd.T_load_out[:, col1] ≈ sr_hyb["load1", :T_load_out]
    end

    # ─────────────────────────────────────────────────────────────────────────
    # :full mode vs :full mode with Vector policy produce the same results
    # ─────────────────────────────────────────────────────────────────────────

    @testset ":full function vs Vector policy consistency" begin
        nw1 = make_network()
        nw2 = make_network()

        # Two policies that produce the same outputs
        mass_flows = 10.0 .+ sin.(2π .* t ./ 3600.0)   # sinusoidal
        temps      = 80.0 .+ 5.0 .* cos.(2π .* t ./ 3600.0)

        pvec = [ProducerOutput(mass_flows[i], temps[i]) for i in 1:N]
        pfn  = (ti, Ta, Tb) -> begin
            i = round(Int, ti / 60) + 1
            ProducerOutput(mass_flows[i], temps[i])
        end

        sr_vec = run_simulation(nw1, t, pvec; mode=:full, ambient_temperature=Tₐ)
        sr_fn  = run_simulation(nw2, t, pfn;  mode=:full, ambient_temperature=Tₐ)

        @test sr_vec.T_load_in      ≈ sr_fn.T_load_in      atol=1e-10
        @test sr_vec.T_producer_out ≈ sr_fn.T_producer_out atol=1e-10
    end

    # ─────────────────────────────────────────────────────────────────────────
    # :full — forward-only recovers same T_load_in as :full
    # (verifying forward step is identical between modes)
    # ─────────────────────────────────────────────────────────────────────────

    @testset ":forward_only T_load_in matches :full T_load_in" begin
        nw1 = make_network()
        nw2 = make_network()

        pvec    = [ProducerOutput(10.0, 80.0) for _ in 1:N]
        p_2arg  = (t, Ta, Tb) -> ProducerOutput(mass_flow=10.0, temperature=80.0)

        sr_full = run_simulation(nw1, t, pvec; mode=:full)
        sr_fwd  = run_simulation(nw2, t, p_2arg; mode=:forward_only)

        # Forward thermal step is identical in both modes
        @test sr_full.T_load_in ≈ sr_fwd.T_load_in atol=1e-10
    end

    # ─────────────────────────────────────────────────────────────────────────
    # Regression: original test preserved
    # ─────────────────────────────────────────────────────────────────────────

    @testset "original regression" begin
        nw = Network()
        j1 = (100.0, 0.0)
        j2 = (250.0, 0.0)
        nw["producer"]  = ProducerNode((0.0, 0.0))
        nw["junction1"] = JunctionNode("J1", j1)
        nw["junction2"] = JunctionNode("J2", j2)
        nw["load1"] = LoadNode("L1", j1 .+ (25*sqrt(2),  25*sqrt(2)), 2.0)
        nw["load2"] = LoadNode("L2", j1 .+ (50*sqrt(2), -50*sqrt(2)), 1.0)
        nw["load3"] = LoadNode("L3", j2 .+ (25*sqrt(2),  25*sqrt(2)), 1.1)
        nw["load4"] = LoadNode("L4", j2 .+ (50*sqrt(2), -50*sqrt(2)), 3.0)
        nw["producer",  "junction1"] = InsulatedPipe(100)
        nw["junction1", "junction2"] = InsulatedPipe(150)
        nw["junction1", "load1"]     = InsulatedPipe(50)
        nw["junction1", "load2"]     = InsulatedPipe(100)
        nw["junction2", "load3"]     = InsulatedPipe(50)
        nw["junction2", "load4"]     = InsulatedPipe(100)

        fill_pipes_with_initial_temperature!(nw, 90.0, 70.0)
        for e in edges(nw)
            @test length(nw[src(e), dst(e)].plugs_f) == 1
            @test nw[src(e), dst(e)].plugs_f[1].T == 90.0
            @test length(nw[src(e), dst(e)].plugs_b) == 1
            @test nw[src(e), dst(e)].plugs_b[1].T == 70.0
        end

        t2 = float.(collect(range(0, stop=3*60*60, step=60)))
        p2 = (t, Tₐ, T_back) -> ProducerOutput(
            mass_flow=15.0, temperature=90 + 10*sin(2π*t/(100*60)))
        results = run_simulation(nw, t2, p2; T0_f=75.0, T0_b=60.0)
        @test length(results) == length(t2)
    end

    # ─────────────────────────────────────────────────────────────────────────
    # SumpNode
    # ─────────────────────────────────────────────────────────────────────────

    @testset "SumpNode in network" begin
        # Topology: producer → (sump) → junction → load1
        #                                         → load2
        # The sump sits on the main supply before the branch.
        nw = Network()
        nw["producer"] = ProducerNode((0.0, 0.0))
        nw["sump"]     = SumpNode("S1", (25.0, 0.0))
        nw["junction"] = JunctionNode((50.0, 0.0))
        nw["load1"]    = LoadNode("L1", (100.0,  50.0), 1.0)
        nw["load2"]    = LoadNode("L2", (100.0, -50.0), 2.0)
        nw["producer", "sump"]     = InsulatedPipe(25)
        nw["sump",     "junction"] = InsulatedPipe(25)
        nw["junction", "load1"]    = InsulatedPipe(50)
        nw["junction", "load2"]    = InsulatedPipe(50)

        # sump_labels should be populated
        @test "sump" ∈ nw.sump_labels
        @test length(nw.sump_labels) == 1

        # replacing a sump with a junction removes it from sump_labels
        nw2 = Network()
        nw2["producer"] = ProducerNode()
        nw2["sump"]     = SumpNode("S")
        nw2["load"]     = LoadNode("L", 1.0)
        nw2["producer", "sump"] = InsulatedPipe(50)
        nw2["sump", "load"]     = InsulatedPipe(50)
        @test "sump" ∈ nw2.sump_labels
        nw2["sump"] = JunctionNode()  # replace with junction
        @test "sump" ∉ nw2.sump_labels

        # removing a sump node removes it from sump_labels
        nw3 = Network()
        nw3["producer"] = ProducerNode()
        nw3["sump"]     = SumpNode()
        nw3["load"]     = LoadNode("L", 1.0)
        nw3["producer", "sump"] = InsulatedPipe(50)
        nw3["sump", "load"]     = InsulatedPipe(50)
        @test "sump" ∈ nw3.sump_labels
        rem_node!(nw3, "sump")
        @test "sump" ∉ nw3.sump_labels
    end

    @testset "SumpNode in SimulationResults" begin
        # Topology: producer → sump → junction → load1 (m_rel=1)
        #                                       → load2 (m_rel=2)
        nw = Network()
        nw["producer"] = ProducerNode((0.0, 0.0))
        nw["sump"]     = SumpNode("S1", (25.0, 0.0))
        nw["junction"] = JunctionNode((50.0, 0.0))
        nw["load1"]    = LoadNode("L1", (100.0,  50.0), 1.0)
        nw["load2"]    = LoadNode("L2", (100.0, -50.0), 2.0)
        nw["producer", "sump"]     = InsulatedPipe(25)
        nw["sump",     "junction"] = InsulatedPipe(25)
        nw["junction", "load1"]    = InsulatedPipe(50)
        nw["junction", "load2"]    = InsulatedPipe(50)

        sim_t  = float.(collect(range(0, stop=60*60, step=60)))
        Nsteps = length(sim_t)
        policy = (t, Ta, Tb) -> ProducerOutput(mass_flow=6.0, temperature=80.0)
        sr = run_simulation(nw, sim_t, policy)

        # sump labels are present
        @test "sump" ∈ sr[:sump_labels]
        @test sr[:sump_labels_dict] isa Dict{String, Int}

        # sump matrices have correct size
        @test size(sr[:mass_flow_sump]) == (Nsteps, 1)
        @test size(sr[:T_sump_f])       == (Nsteps, 1)
        @test size(sr[:T_sump_b])       == (Nsteps, 1)

        # convenience indexing
        mf   = sr["sump", :mass_flow_sump]
        T_f  = sr["sump", :T_sump_f]
        T_b  = sr["sump", :T_sump_b]
        @test length(mf)  == Nsteps
        @test length(T_f) == Nsteps
        @test length(T_b) == Nsteps

        # mass flow at sump must equal total producer mass flow (single path before branch)
        @test all(isapprox.(mf, 6.0; atol=1e-10))

        # forward temperature should be bounded (between initial T and supply T)
        @test all(T_f[2:end] .>= 25.0)   # above initial pipe temperature
        @test all(T_f[2:end] .<= 80.0)   # at most the producer supply temperature

        # return temperature should be lower than supply
        @test all(T_b[2:end] .>= 25.0)
        @test all(T_b[2:end] .<= 80.0)
    end

    @testset "SumpNode in forward_only mode" begin
        nw = Network()
        nw["producer"] = ProducerNode()
        nw["sump"]     = SumpNode()
        nw["load"]     = LoadNode("L", 1.0)
        nw["producer", "sump"] = InsulatedPipe(50)
        nw["sump",     "load"] = InsulatedPipe(50)

        sim_t = float.(collect(0:60.0:3600.0))
        N     = length(sim_t)
        policy = (t, Ta, Tb) -> ProducerOutput(mass_flow=5.0, temperature=75.0)
        sr = run_simulation(nw, sim_t, policy; mode=:forward_only)

        # T_sump_b is Nothing in forward_only mode
        @test isnothing(sr[:T_sump_b])
        @test isnothing(sr["sump", :T_sump_b])

        # T_sump_f is recorded
        @test size(sr[:T_sump_f]) == (N, 1)
        @test !all(isnan, sr["sump", :T_sump_f])
    end

    @testset "SimulationResults show and length with sumps" begin
        nw = make_network()  # no sumps
        sim_t  = float.(collect(0:60.0:3600.0))
        policy = (t, Ta, Tb) -> ProducerOutput(mass_flow=6.0, temperature=80.0)
        sr = run_simulation(nw, sim_t, policy)

        @test length(sr.sump_labels) == 0
        @test size(sr[:mass_flow_sump]) == (length(sim_t), 0)
        @test contains(string(sr), "0 sump node(s)")
    end

    # ─────────────────────────────────────────────────────────────────────────
    # Time-varying m_rel
    # ─────────────────────────────────────────────────────────────────────────

    @testset "time-varying m_rel: validation errors" begin
        # mixing constant and vector m_rel across loads is not allowed
        nw_mix = make_network()
        set_load_m_rel!(nw_mix, "load1", fill(1.0, N))   # vector
        # load2 still has scalar m_rel=2.0 → mixed → must error
        @test_throws ErrorException run_simulation(nw_mix, t,
            (t, Ta, Tb) -> ProducerOutput(mass_flow=10.0, temperature=80.0))

        # vector m_rel with wrong length
        nw_bad = make_network()
        set_load_m_rel!(nw_bad, "load1", fill(1.0, N - 1))   # wrong length
        set_load_m_rel!(nw_bad, "load2", fill(2.0, N - 1))
        @test_throws ErrorException run_simulation(nw_bad, t,
            (t, Ta, Tb) -> ProducerOutput(mass_flow=10.0, temperature=80.0))
    end

    @testset "time-varying m_rel: constant result when vector is uniform" begin
        # A time-varying m_rel vector that is constant should give the same
        # simulation results as the equivalent scalar m_rel.
        nw_const  = make_network()  # load1 m_rel=1, load2 m_rel=2
        nw_vector = make_network()
        set_load_m_rel!(nw_vector, "load1", fill(1.0, N))
        set_load_m_rel!(nw_vector, "load2", fill(2.0, N))

        policy = (ts, Ta, Tb) -> ProducerOutput(mass_flow=10.0, temperature=80.0)
        sr_c = run_simulation(nw_const,  t, policy; mode=:forward_only)
        sr_v = run_simulation(nw_vector, t, policy; mode=:forward_only)

        @test sr_c.mass_flow_load ≈ sr_v.mass_flow_load   atol=1e-10
        @test sr_c.T_load_in      ≈ sr_v.T_load_in        atol=1e-10
    end

    @testset "time-varying m_rel: flows follow the prescribed vector" begin
        # Two-step pattern: odd steps → equal split (1:1), even steps → 3:1 split.
        # Verify that the recorded mass flows at each step match expectations.
        total_flow = 12.0
        m1 = [isodd(i) ? 1.0 : 3.0 for i in 1:N]
        m2 = [isodd(i) ? 1.0 : 1.0 for i in 1:N]

        nw_tv = make_network()
        set_load_m_rel!(nw_tv, "load1", m1)
        set_load_m_rel!(nw_tv, "load2", m2)

        policy = (ts, Ta, Tb) -> ProducerOutput(mass_flow=total_flow, temperature=80.0)
        sr = run_simulation(nw_tv, t, policy; mode=:forward_only)

        # Check the first two steps (indices 1 and 2) against expected fractions.
        # step 1: m_rel = (1, 1) → each load gets 0.5 of total
        expected_step1 = total_flow * 1.0 / (1.0 + 1.0)
        @test sr.mass_flow_load[1, sr[:load_labels_dict]["load1"]] ≈ expected_step1  atol=1e-10
        @test sr.mass_flow_load[1, sr[:load_labels_dict]["load2"]] ≈ expected_step1  atol=1e-10

        # step 2: m_rel = (3, 1) → load1 gets 3/4, load2 gets 1/4
        expected_step2_l1 = total_flow * 3.0 / (3.0 + 1.0)
        expected_step2_l2 = total_flow * 1.0 / (3.0 + 1.0)
        @test sr.mass_flow_load[2, sr[:load_labels_dict]["load1"]] ≈ expected_step2_l1  atol=1e-10
        @test sr.mass_flow_load[2, sr[:load_labels_dict]["load2"]] ≈ expected_step2_l2  atol=1e-10
    end

    @testset "time-varying m_rel: full mode completes without error" begin
        nw_tv = make_network()
        set_load_m_rel!(nw_tv, "load1", fill(1.0, N))
        set_load_m_rel!(nw_tv, "load2", fill(2.0, N))
        policy = (ts, Ta, Tb) -> ProducerOutput(mass_flow=10.0, temperature=80.0)
        sr = run_simulation(nw_tv, t, policy; mode=:full, ambient_temperature=Tₐ)
        @test length(sr) == N
        @test !any(isnan, sr.T_load_in)
        @test !any(isnan, sr.T_load_out)
    end

end
