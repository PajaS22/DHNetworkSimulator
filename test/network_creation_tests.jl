@testset "Network creation helpers" begin

    @testset "Network from DiGraph" begin
        g = SimpleDiGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 2, 4)
        nw = Network(g)
        @test nv(nw) == 4
        @test ne(nw) == 3
        # nodes and edges are placeholder types
        @test all(n isa EmptyNode for n in vertices_data(nw))
        @test all(e isa EmptyEdge for e in edges_data(nw))
        # default labels are "1", "2", ...
        @test has_label(nw, "1")
        @test has_label(nw, "4")
    end

    @testset "Network from DiGraph rejects invalid graphs" begin
        # cyclic graph
        g_cyclic = SimpleDiGraph(3)
        add_edge!(g_cyclic, 1, 2)
        add_edge!(g_cyclic, 2, 3)
        add_edge!(g_cyclic, 3, 1)
        @test_throws ErrorException Network(g_cyclic)

        # disconnected graph
        g_disconnected = SimpleDiGraph(3)
        add_edge!(g_disconnected, 1, 2)
        # node 3 is isolated — not connected to 1 or 2
        @test_throws ErrorException Network(g_disconnected)
    end

    @testset "name_nodes!" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        nw = Network(g)
        name_nodes!(nw, ["root", "leaf1", "leaf2"])
        @test has_label(nw, "root")
        @test has_label(nw, "leaf1")
        @test has_label(nw, "leaf2")
        @test !has_label(nw, "1")
        @test !has_label(nw, "2")
        @test !has_label(nw, "3")
        @test nv(nw) == 3
        @test ne(nw) == 2

        # wrong number of labels → error
        g2 = SimpleDiGraph(2)
        add_edge!(g2, 1, 2)
        nw2 = Network(g2)
        @test_throws ErrorException name_nodes!(nw2, ["only_one"])
    end

    @testset "identify_producer_and_loads!" begin
        g = SimpleDiGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 2, 4)
        nw = Network(g)
        name_nodes!(nw, ["producer", "junction", "load1", "load2"])
        identify_producer_and_loads!(nw)

        @test nw["producer"] isa ProducerNode
        @test nw["junction"] isa JunctionNode
        @test nw["load1"] isa LoadNode
        @test nw["load2"] isa LoadNode
        @test nw.producer_label == "producer"
        @test nw.load_labels == Set(["load1", "load2"])
    end

    @testset "fill_physical_params!" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        nw = Network(g)
        params = Dict(
            (1, 2) => PipeParams(; length=120.0, inner_diameter=0.08),
            (1, 3) => PipeParams(; length=80.0,  inner_diameter=0.06),
        )
        fill_physical_params!(nw, params)

        @test nw["1", "2"] isa InsulatedPipe
        @test nw["1", "3"] isa InsulatedPipe
        @test length(nw["1", "2"]) == 120.0
        @test inner_diameter(nw["1", "2"]) == 0.08
        @test length(nw["1", "3"]) == 80.0
    end

    @testset "fill_load_specs! per-node" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        nw = Network(g)
        name_nodes!(nw, ["src", "l1", "l2"])
        identify_producer_and_loads!(nw)

        # coefficients chosen so polynomial_load is ≥ 0 over [-30, 30]
        pwr = Dict("l1" => (500.0, -30.0, 0.5), "l2" => (300.0, -10.0, 0.3))
        mr  = Dict("l1" => 1.0, "l2" => 2.0)
        fill_load_specs!(nw, pwr, mr)

        @test nw["l1"].load isa LoadSpec
        @test nw["l1"].load.fn === polynomial_load
        @test nw["l1"].load.params == [500.0, -30.0, 0.5]
        @test nw["l1"].m_rel == 1.0
        @test nw["l2"].load.params == [300.0, -10.0, 0.3]
        @test nw["l2"].m_rel == 2.0
    end

    @testset "identify_sumps!" begin
        # convert EmptyNode to SumpNode
        nw = Network()
        nw["producer"] = ProducerNode()
        nw["mid"]      = EmptyNode()
        nw["load"]     = LoadNode("L", 1.0)
        nw["producer", "mid"]  = InsulatedPipe(50)
        nw["mid",      "load"] = InsulatedPipe(50)

        identify_sumps!(nw, ["mid"])
        @test nw["mid"] isa SumpNode
        @test "mid" ∈ nw.sump_labels
        @test nw["mid"].common.info == "mid"

        # convert JunctionNode to SumpNode — info and position are preserved
        nw2 = Network()
        nw2["producer"] = ProducerNode()
        nw2["jct"]      = JunctionNode("my junction", (5.0, 6.0))
        nw2["load1"]    = LoadNode("L1", 1.0)
        nw2["load2"]    = LoadNode("L2", 1.0)
        nw2["producer", "jct"]   = InsulatedPipe(50)
        nw2["jct",      "load1"] = InsulatedPipe(50)
        nw2["jct",      "load2"] = InsulatedPipe(50)

        identify_sumps!(nw2, ["jct"])
        @test nw2["jct"] isa SumpNode
        @test "jct" ∈ nw2.sump_labels
        @test nw2["jct"].common.info     == "my junction"   # info preserved
        @test nw2["jct"].common.position == (5.0, 6.0)      # position preserved

        # already a SumpNode → silently skipped
        identify_sumps!(nw2, ["jct"])
        @test nw2["jct"] isa SumpNode

        # single-label convenience form
        nw3 = Network()
        nw3["producer"] = ProducerNode()
        nw3["s"]        = JunctionNode("S")
        nw3["load"]     = LoadNode("L", 1.0)
        nw3["producer", "s"]    = InsulatedPipe(50)
        nw3["s",        "load"] = InsulatedPipe(50)
        identify_sumps!(nw3, "s")
        @test nw3["s"] isa SumpNode

        # error on non-existent label
        @test_throws ErrorException identify_sumps!(nw, ["nonexistent"])

        # error on incompatible node types
        nw4 = Network()
        nw4["producer"] = ProducerNode()
        nw4["load"]     = LoadNode("L", 1.0)
        nw4["producer", "load"] = InsulatedPipe(50)
        @test_throws ErrorException identify_sumps!(nw4, ["producer"])
        @test_throws ErrorException identify_sumps!(nw4, ["load"])
    end

    @testset "fill_load_specs! default" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        nw = Network(g)
        name_nodes!(nw, ["src", "la", "lb"])
        identify_producer_and_loads!(nw)

        fill_load_specs!(nw; pwr_coefs=(200.0, 0.0, 0.0), m_r=1.5)

        @test nw["la"].load isa LoadSpec
        @test nw["la"].load.params == [200.0, 0.0, 0.0]
        @test nw["la"].m_rel == 1.5
        @test nw["lb"].load.params == [200.0, 0.0, 0.0]
        @test nw["lb"].m_rel == 1.5
    end

end

# helper to build a small two-load network for LoadSpec API tests
function _two_load_network()
    g = SimpleDiGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 1, 3)
    nw = Network(g)
    name_nodes!(nw, ["src", "la", "lb"])
    identify_producer_and_loads!(nw)
    return nw
end

@testset "LoadSpec API" begin

    @testset "validate_load_spec — accepts valid function" begin
        @test_nowarn validate_load_spec(polynomial_load, [540.0, -36.0, 0.6])
        @test_nowarn validate_load_spec(polynomial_load, [200.0, 0.0, 0.0])  # constant
        custom_fn = (p, T) -> p[1] * exp(-p[2] * (T - p[3])^2)  # Gaussian, always ≥ 0
        @test_nowarn validate_load_spec(custom_fn, [500.0, 0.01, -10.0])
    end

    @testset "validate_load_spec — rejects negative values" begin
        # linear with negative slope goes negative before T_a = 30
        @test_throws ErrorException validate_load_spec(polynomial_load, [100.0, -5.0, 0.0])
    end

    @testset "validate_load_spec — rejects non-finite values" begin
        bad_fn = (p, T) -> T == 0.0 ? Inf : p[1]
        @test_throws ErrorException validate_load_spec(bad_fn, [100.0])
    end

    @testset "set_load_fn! — single load" begin
        nw = _two_load_network()
        set_load_fn!(nw, "la", polynomial_load, [400.0, -20.0, 0.4])
        @test nw["la"].load isa LoadSpec
        @test nw["la"].load.fn === polynomial_load
        @test nw["la"].load.params == [400.0, -20.0, 0.4]
    end

    @testset "set_load_fn! — rejects invalid function" begin
        nw = _two_load_network()
        @test_throws ErrorException set_load_fn!(nw, "la", polynomial_load, [100.0, -5.0, 0.0])
    end

    @testset "set_load_fn! — non-existent or wrong node type" begin
        nw = _two_load_network()
        @test_throws ErrorException set_load_fn!(nw, "does_not_exist", polynomial_load, [100.0])
        @test_throws ErrorException set_load_fn!(nw, "src", polynomial_load, [100.0])  # producer, not load
    end

    @testset "set_load_fn! — all loads, uniform params" begin
        nw = _two_load_network()
        set_load_fn!(nw, polynomial_load, [540.0, -36.0, 0.6])
        @test nw["la"].load.params == [540.0, -36.0, 0.6]
        @test nw["lb"].load.params == [540.0, -36.0, 0.6]
        # nodes hold independent copies
        nw["la"].load.params[1] = 999.0
        @test nw["lb"].load.params[1] == 540.0
    end

    @testset "set_load_fn! — all loads, per-load params dict" begin
        nw = _two_load_network()
        custom_fn = (p, T) -> p[1]  # constant, always ≥ 0
        params_dict = Dict("la" => [300.0], "lb" => [500.0])
        set_load_fn!(nw, custom_fn, params_dict)
        @test nw["la"].load.fn === custom_fn
        @test nw["la"].load.params == [300.0]
        @test nw["lb"].load.fn === custom_fn
        @test nw["lb"].load.params == [500.0]
    end

    @testset "set_load_fn! dict — atomic: no partial update on error" begin
        nw = _two_load_network()
        set_load_fn!(nw, "la", polynomial_load, [540.0, -36.0, 0.6])  # la has a valid spec
        bad_dict = Dict("la" => [999.0], "lb" => [100.0, -5.0, 0.0])  # lb params invalid
        @test_throws ErrorException set_load_fn!(nw, polynomial_load, bad_dict)
        # la must be unchanged since validation runs before any update
        @test nw["la"].load.params == [540.0, -36.0, 0.6]
    end

    @testset "set_load_params! (Vector) — updates params, keeps function" begin
        nw = _two_load_network()
        set_load_fn!(nw, "la", polynomial_load, [540.0, -36.0, 0.6])
        set_load_params!(nw, "la", [400.0, -20.0, 0.4])
        @test nw["la"].load.fn === polynomial_load   # function unchanged
        @test nw["la"].load.params == [400.0, -20.0, 0.4]
    end

    @testset "set_load_params! (Vector) — rejects invalid params" begin
        nw = _two_load_network()
        set_load_fn!(nw, "la", polynomial_load, [540.0, -36.0, 0.6])
        @test_throws ErrorException set_load_params!(nw, "la", [100.0, -5.0, 0.0])
    end

    @testset "set_load_params! (Vector) — errors when no function set" begin
        nw = _two_load_network()  # loads have default LoadSpec from constructor
        # default spec is valid; explicitly set to missing to test guard
        nw["la"].load = missing
        @test_throws ErrorException set_load_params!(nw, "la", [540.0, -36.0, 0.6])
    end

    @testset "set_load_params! (Dict) — bulk update" begin
        nw = _two_load_network()
        set_load_fn!(nw, polynomial_load, [540.0, -36.0, 0.6])
        new_params = Dict("la" => [400.0, -20.0, 0.4], "lb" => [300.0, -15.0, 0.3])
        set_load_params!(nw, new_params)
        @test nw["la"].load.params == [400.0, -20.0, 0.4]
        @test nw["lb"].load.params == [300.0, -15.0, 0.3]
    end

    @testset "power_consumption uses LoadSpec" begin
        node = LoadNode("test"; load=LoadSpec(polynomial_load, [540.0, -36.0, 0.6]))
        # at T_a = 0: P = 540 kW → 540_000 W
        @test power_consumption(node, 0.0) ≈ 540_000.0
        # at T_a = 30: P ≈ 0 kW → clamped to 0 W
        @test power_consumption(node, 30.0) ≥ 0.0
        # custom constant function
        node2 = LoadNode("test2"; load=LoadSpec((p, T) -> p[1], [200.0]))
        @test power_consumption(node2, -15.0) ≈ 200_000.0
        @test power_consumption(node2, 25.0)  ≈ 200_000.0
    end

end
