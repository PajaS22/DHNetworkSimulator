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

        pwr = Dict("l1" => (500.0, -30.0, 0.5), "l2" => (300.0, -20.0, 0.3))
        mr  = Dict("l1" => 1.0, "l2" => 2.0)
        fill_load_specs!(nw, pwr, mr)

        @test nw["l1"].load == (500.0, -30.0, 0.5)
        @test nw["l1"].m_rel == 1.0
        @test nw["l2"].load == (300.0, -20.0, 0.3)
        @test nw["l2"].m_rel == 2.0
    end

    @testset "fill_load_specs! default" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        nw = Network(g)
        name_nodes!(nw, ["src", "la", "lb"])
        identify_producer_and_loads!(nw)

        fill_load_specs!(nw; pwr_coefs=(100.0, -5.0, 0.0), m_r=1.5)

        @test nw["la"].load == (100.0, -5.0, 0.0)
        @test nw["la"].m_rel == 1.5
        @test nw["lb"].load == (100.0, -5.0, 0.0)
        @test nw["lb"].m_rel == 1.5
    end

end
