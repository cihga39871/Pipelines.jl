
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

    @test_nowarn p = @pkg JuliaProgram(
        id_file = "id_file",
        inputs = ["a",
                  "b" => Int],
        outputs = "c" => "<a>.<b>",
        main = quote
            println("inputs are ", a, " and ", b)
            println("You can also use info in outputs: ", c)
            println("The returned value will be assigned to a new outputs")
            c = b^2
        end)
end
