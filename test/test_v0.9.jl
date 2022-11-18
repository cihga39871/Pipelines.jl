
@testset "v0.9 rerun" begin
    int_string(file) = parse(Int, chomp(String(read(file))))
    jp = JuliaProgram(
        name = "Test v0.9 Rerun",
        id_file = ".v09rerun",
        inputs = [
            "infile",
            "add" => Int
        ],
        outputs = [
            "outfile" => "<infile>.out"
        ],
        main = quote
            x = $int_string(infile)
            open(outfile, "w+") do io
                print(io, x + add)
            end
        end
    )
    
    infile = "in.txt"
    open(infile, "w+") do io
        print(io, 10)
    end
    run(jp, infile = infile, add = 3)

    @test_warn "Skipped finished program" run(jp, infile = infile, add = 3)
    outfile = "in.txt.out"
    @test int_string(outfile) == 13

    run(jp, infile = infile, add = 5)
    @test int_string(outfile) == 15

    run(jp, infile = infile, add = 3)
    @test int_string(outfile) == 13

    open(infile, "w+") do io
        print(io, 100)
    end
    run(jp, infile = infile, add = 3)
    @test int_string(outfile) == 103
end

@testset "v0.9 quote function dictreplace!" begin
    inputs = ["abc", "def"]
    f(abc; def = "123") = @show abc def
    expr = quote
        $f(abc, def = def)
    end
    f2 = quote_function(expr, inputs)
    f2(Dict("abc" => 2, "def" => 33))
end