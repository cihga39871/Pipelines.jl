
@testset "Arg" begin

    @test Arg("argname", 5) == Arg("argname", Any, 5, false, false)
    @test Arg("argname", Real) == Arg("argname", Real, nothing, true, false)
    @test Arg(:argname, Real, 2.5) == Arg("argname", Real, 2.5, false, true)
    @test_throws ErrorException Arg("argname", Int, 0.6)
    @test_throws ErrorException Arg("name", Int, 0.6)

    @test Arg("x") == Arg("x", Any, nothing)
    @test Arg("x" => Int) == Arg("x", Int, nothing)
    @test Arg("x" => 5) == Arg("x", Any, 5)
    @test Arg("x" => 5 => Int) == Arg("x", Int, 5)
    @test Arg("x" => Int => 5) == Arg("x", Int, 5)
    @test Arg("x" => Int => DataType) == Arg("x", DataType, Int)
    @test Arg("x" => DataType => Int) == Arg("x", DataType, Int)
    @test Arg( :b => (5.0 => Int)) == Arg("b", Int, 5, false, true)

    @test_throws ErrorException Arg("x" => Int => "5")
    @test_throws ErrorException Arg("x" => "5" => Int)
    @test_throws ErrorException Arg(5)

    v = ["A", "B" => Int, "C" => 9.0 => Float64]
    @test Pipelines.parse_arg(v) == Arg[Arg(v[1]), Arg(v[2]), Arg(v[3])]

    v = ["A", :B => Int, "C" => 9.0 => Float64]
    @test Pipelines.parse_arg(v) == Arg[Arg(v[1]), Arg(v[2]), Arg(v[3])]

    v = ["A", :B => Int => 5.0, "C" => 9.0 => Float64]
    @test Pipelines.parse_arg(v) == Arg[Arg(v[1]), Arg(v[2]), Arg(v[3])]

    @test_nowarn display(Arg[])

end

@testset "Arg forward" begin
    @test_throws "cannot use the reserved name" JuliaProgram(
        id_file = "id_file",
        inputs = [
            "a",
            "b" => Int,
            "NAME", :thread, :ncpu
        ],
        outputs = [
            "c" => "<a>.<b>"
        ],
        main = quote
            c = "$a.$b.$NAME.$thread"
        end,
        arg_forward = "NAME" => "name"
    )
    p = JuliaProgram(
        id_file = "id_file",
        inputs = [
            "a",
            "b" => Int,
            "NAME", :thread
        ],
        outputs = [
            "c" => "<a>.<b>"
        ],
        main = quote
            c = "$a.$b.$NAME.$thread"
        end,
        arg_forward = "NAME" => "name"
    )
    @test p.arg_forward == ["NAME" => :name]

    p = JuliaProgram(
        id_file = "id_file",
        inputs = [
            "a",
            "b" => Int,
            "NAME", :thread
        ],
        outputs = [
            "c" => "<a>.<b>"
        ],
        main = quote
            c = "$a.$b.$NAME.$thread"
        end,
        arg_forward = ["NAME" => "name", :thread => :ncpu]
    )
    @test p.arg_forward == ["NAME" => :name, "thread" => :ncpu]
    p = JuliaProgram(
        id_file = "id_file",
        inputs = [
            "a",
            "b" => Int,
            "NAME", :thread
        ],
        outputs = [
            "c" => "<a>.<b>"
        ],
        main = quote
            c = "$a.$b.$NAME.$thread"
        end,
        arg_forward = Dict("NAME" => "name", :thread => :ncpu)
    )
    @test p.arg_forward == ["NAME" => :name, "thread" => :ncpu] || 
    p.arg_forward == ["thread" => :ncpu, "NAME" => :name]

end