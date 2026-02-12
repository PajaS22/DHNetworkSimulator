@testset "Network" begin
    @testset "creation" begin
        nw = Network()
        @test nv(nw) == 0
        @test ne(nw) == 0
        @test isnothing(nw.producer_label)
        @test isempty(nw.load_labels)
    end
    @testset "adding nodes and edges" begin
        nw = Network()
        nw["node1"] = EmptyNode()
        nw["node2"] = EmptyNode()
        @test nv(nw) == 2

        nw["node1", "node2"] = EmptyEdge()
        
        @test ne(nw) == 1
    end
    @testset "indexing" begin
        nw = Network()
        nw["node1"] = EmptyNode()
        nw["node2"] = EmptyNode()
        nw["node1", "node2"] = EmptyEdge()

        n1 = index_for(nw, "node1")
        n2 = index_for(nw, "node2")
        @test n1 != n2

        @test nw[n1] === nw["node1"]
        @test nw[n2] === nw["node2"]

        e = nw[n1, n2]
        @test e === nw["node1", "node2"]
    end
    @testset "removing nodes" begin
        nw = Network()
        nw["node1"] = EmptyNode()
        nw["node2"] = EmptyNode()
        nw["node1", "node2"] = EmptyEdge()

        rem_node!(nw, "node1")
        @test nv(nw) == 1
        @test ne(nw) == 0
        @test !has_label(nw, "node1")

        nw1 = Network()
        nw1["node1"] = EmptyNode()
        nw1["node2"] = EmptyNode()
        rem_node!(nw1, index_for(nw1, "node2"))
        @test nv(nw1) == 1
        @test !has_label(nw1, "node2")
    end
    @testset "vertices_data and edges_data" begin
        nw = Network()
        nw["node1"] = EmptyNode()
        nw["node2"] = EmptyNode()
        nw["node1", "node2"] = EmptyEdge()

        vdata = vertices_data(nw)
        @test length(vdata) == 2

        edata = edges_data(nw)
        @test length(edata) == 1
    end
    @testset "neighbors and degree functions" begin
        nw = Network()
        nw["node1"] = EmptyNode()
        nw["node2"] = EmptyNode()
        nw["node3"] = EmptyNode()
        nw["node1", "node2"] = EmptyEdge()
        nw["node1", "node3"] = EmptyEdge()

        n1 = index_for(nw, "node1")
        n2 = index_for(nw, "node2")
        n3 = index_for(nw, "node3")

        l1 = label_for(nw, n1)
        l2 = label_for(nw, n2)
        l3 = label_for(nw, n3)

        @test DHN.degree(nw, "node1") == 2
        @test DHN.outdegree(nw, "node1") == 2
        @test DHN.indegree(nw, "node1") == 0
        @test DHN.inneighbors(nw, "node2") == ["node1"]
    end
    @testset "load labels" begin
        # keep track of load node ids when adding, removing and replacing load nodes
        nw = Network()
        nw["load1"] = LoadNode()
        nw["load2"] = LoadNode()
        nw["load3"] = LoadNode()
        nw["junction"] = JunctionNode()

        @test nw.load_labels == Set(["load1", "load2", "load3"])

        # remove load1
        rem_node!(nw, "load1")
        @test nw.load_labels == Set(["load2", "load3"])
        # replace load2 with a new load node
        nw["load2"] = LoadNode()
        @test nw.load_labels == Set(["load2", "load3"])

        # replace load2 with a junction node
        nw["load2"] = JunctionNode()
        @test nw.load_labels == Set(["load3"])
        # replace load3 with a producer node
        nw["load3"] = ProducerNode()
        @test isempty(nw.load_labels)

        # add a new load node
        nw["load4"] = LoadNode()
        @test nw.load_labels == Set(["load4"])
        # replace load4 with an Empty node
        nw["load4"] = EmptyNode()
        @test isempty(nw.load_labels)
    end
    @testset "producer label" begin
        # keep track of producer node id when adding, removing and replacing producer node
        nw = Network()
        @test isnothing(nw.producer_label)

        nw["producer"] = ProducerNode()
        @test nw.producer_label == "producer"

        # replace producer with a new producer node
        nw["producer"] = ProducerNode()
        @test nw.producer_label == "producer"
        # replace producer with a junction node
        nw["producer"] = JunctionNode()
        @test isnothing(nw.producer_label)

        # add a new producer node
        nw["new_producer"] = ProducerNode()
        @test nw.producer_label == "new_producer"

        # remove producer node
        rem_node!(nw, "new_producer")
        @test isnothing(nw.producer_label)
        
        nw["producer"] = ProducerNode() # change Junction to Producer
        @test nw.producer_label == "producer"

        # replace producer with a load node
        nw["producer"] = LoadNode()
        @test isnothing(nw.producer_label)

        # replace producer with an Empty node
        nw["producer"] = ProducerNode() # change Load to Producer
        nw["producer"] = EmptyNode() # change Producer to Empty
        @test isnothing(nw.producer_label)

        # try adding second producer node (should error)
        nw["producer1"] = ProducerNode()
        @test nw.producer_label == "producer1"

        try
            nw["producer2"] = ProducerNode()
            @test false  # Pokud k chybě nedojde, test selže
        catch e
            @test e isa ErrorException
        end
    end


end