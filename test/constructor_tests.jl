@testset "Constructors" begin

    @testset "PipeParams" begin
        # keyword constructor with defaults
        p = PipeParams(; length=100.0, inner_diameter=0.1)
        @test p.length == 100.0
        @test p.inner_diameter == 0.1
        @test p.heat_resistance_forward == 3.0
        @test p.heat_resistance_backward == 4.0

        # keyword constructor with all args
        p2 = PipeParams(; length=50.0, inner_diameter=0.2,
                          heat_resistance_forward=2.0, heat_resistance_backward=3.0)
        @test p2.heat_resistance_forward == 2.0
        @test p2.heat_resistance_backward == 3.0

        # Int inputs are promoted to Float64
        p3 = PipeParams(; length=80, inner_diameter=1)
        @test p3.length === 80.0
        @test p3.inner_diameter === 1.0
    end

    @testset "InsulatedPipe" begin
        # all-keyword constructor with defaults
        pipe = InsulatedPipe()
        @test pipe.info == "pipe"
        @test length(pipe) == 100.0
        @test pipe_length(pipe) == length(pipe)  # pipe_length is an alias for length
        @test inner_diameter(pipe) == 0.1
        @test heat_resistance_forward(pipe) == 3.0
        @test heat_resistance_backward(pipe) == 4.0
        @test ismissing(pipe.mass_flow)
        @test ismissing(pipe.m_rel)
        @test isempty(pipe.plugs_f)
        @test isempty(pipe.plugs_b)

        # all-keyword with explicit values
        pipe2 = InsulatedPipe(; info="main", length=200.0, inner_diameter=0.15,
                                heat_resistance_forward=2.5, heat_resistance_backward=3.5,
                                mass_flow=5.0, m_rel=0.6)
        @test pipe2.info == "main"
        @test length(pipe2) == 200.0
        @test inner_diameter(pipe2) == 0.15
        @test heat_resistance_forward(pipe2) == 2.5
        @test heat_resistance_backward(pipe2) == 3.5
        @test pipe2.mass_flow == 5.0
        @test pipe2.m_rel == 0.6

        # Int length is promoted to Float64
        pipe3 = InsulatedPipe(; length=150)
        @test length(pipe3) === 150.0

        # from PipeParams (default info)
        params = PipeParams(; length=80.0, inner_diameter=0.08)
        pipe4 = InsulatedPipe(params)
        @test pipe4.info == "pipe"
        @test length(pipe4) == 80.0
        @test inner_diameter(pipe4) == 0.08

        # from PipeParams with m_rel keyword
        pipe5 = InsulatedPipe(params; m_rel=0.7)
        @test pipe5.m_rel == 0.7

        # info + PipeParams
        pipe6 = InsulatedPipe("branch", params)
        @test pipe6.info == "branch"
        @test length(pipe6) == 80.0

        # info + PipeParams with keywords
        pipe7 = InsulatedPipe("supply", params; mass_flow=10.0, m_rel=1.0)
        @test pipe7.mass_flow == 10.0
        @test pipe7.m_rel == 1.0

        # length shorthand (Int)
        pipe8 = InsulatedPipe(300)
        @test length(pipe8) === 300.0
        @test pipe8.info == "pipe"

        # length shorthand (Float64)
        pipe9 = InsulatedPipe(75.5)
        @test length(pipe9) === 75.5
    end

    @testset "ZeroPipe" begin
        z = ZeroPipe()
        @test length(z) == 0.0
        @test inner_diameter(z) == 0.0
        @test heat_resistance_forward(z) == 0.0
        @test heat_resistance_backward(z) == 0.0
        @test z.info == "zero pipe"
        @test ismissing(z.mass_flow)
        @test ismissing(z.m_rel)

        z2 = ZeroPipe("connector")
        @test z2.info == "connector"
        @test length(z2) == 0.0

        z3 = ZeroPipe(; mass_flow=2.5, m_rel=1.0)
        @test z3.mass_flow == 2.5
        @test z3.m_rel == 1.0

        z4 = ZeroPipe("named"; m_rel=0.5)
        @test z4.info == "named"
        @test z4.m_rel == 0.5
    end

    @testset "ZeroPipe identity" begin
        @test ZeroPipe() isa ZeroPipe
        @test ZeroPipe("x") isa ZeroPipe
        @test !(InsulatedPipe(100) isa ZeroPipe)
    end

    @testset "JunctionNode constructors" begin
        j1 = JunctionNode()
        @test j1.common.info == "junction"
        @test ismissing(j1.common.position)

        j2 = JunctionNode("J2")
        @test j2.common.info == "J2"

        j3 = JunctionNode((3.0, 4.0))
        @test j3.common.info == "junction"
        @test j3.common.position == (3.0, 4.0)

        j4 = JunctionNode("main junction", (1.0, 2.0))
        @test j4.common.info == "main junction"
        @test j4.common.position == (1.0, 2.0)
    end

    @testset "ProducerNode constructors" begin
        p1 = ProducerNode()
        @test p1.common.info == "producer"
        @test ismissing(p1.common.position)

        p2 = ProducerNode("source")
        @test p2.common.info == "source"

        p3 = ProducerNode((0.0, 0.0))
        @test p3.common.position == (0.0, 0.0)

        p4 = ProducerNode("P", (1.0, 2.0))
        @test p4.common.info == "P"
        @test p4.common.position == (1.0, 2.0)
    end

    @testset "SumpNode constructors" begin
        s1 = SumpNode()
        @test s1.common.info == "sump"
        @test ismissing(s1.common.position)

        s2 = SumpNode("S2")
        @test s2.common.info == "S2"

        s3 = SumpNode((7.0, 8.0))
        @test s3.common.info == "sump"
        @test s3.common.position == (7.0, 8.0)

        s4 = SumpNode("main sump", (1.0, 2.0))
        @test s4.common.info == "main sump"
        @test s4.common.position == (1.0, 2.0)

        @test SumpNode() isa NodeType
        @test SumpNode() isa SumpNode
    end

    @testset "LoadNode constructors" begin
        l1 = LoadNode()
        @test l1.common.info == "load"
        @test ismissing(l1.common.position)
        @test ismissing(l1.m_rel)

        l2 = LoadNode("consumer")
        @test l2.common.info == "consumer"

        l3 = LoadNode((5.0, 5.0))
        @test l3.common.info == "load"
        @test l3.common.position == (5.0, 5.0)

        l4 = LoadNode("L4", (1.0, 2.0))
        @test l4.common.info == "L4"
        @test l4.common.position == (1.0, 2.0)
        @test ismissing(l4.m_rel)

        # with m_rel as 3rd positional arg
        l5 = LoadNode("L5", (0.0, 0.0), 0.5)
        @test l5.m_rel == 0.5

        # with custom LoadSpec as 3rd positional arg
        custom_spec = LoadSpec(polynomial_load, [100.0, -5.0, 0.1])
        l6 = LoadNode("L6", (0.0, 0.0), custom_spec)
        @test l6.load.fn === polynomial_load
        @test l6.load.params == [100.0, -5.0, 0.1]
        @test ismissing(l6.m_rel)

        # with keyword load override
        custom_spec2 = LoadSpec(polynomial_load, [200.0, -10.0, 0.2])
        l7 = LoadNode("L7"; load=custom_spec2)
        @test l7.load.params == [200.0, -10.0, 0.2]

        # default load spec uses polynomial_load and DEFAULT_LOAD_PARAMS
        @test l1.load isa LoadSpec
        @test l1.load.fn === polynomial_load
        @test l1.load.params == DEFAULT_LOAD_PARAMS

        # each node gets an independent copy of the default params
        l1.load.params[1] = 999.0
        l8 = LoadNode()
        @test l8.load.params[1] == DEFAULT_LOAD_PARAMS[1]  # l8 unaffected
    end

end
