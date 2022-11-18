
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

@testset "v0.9.3 Possible files in Cmd also write to run_id_file" begin

    prog_cmd = JuliaProgram(
        name = "Program Cmd",
        id_file = ".cmd",
        inputs = "CMD" => Base.AbstractCmd,
        main = quote
            run(CMD)
        end
    )

    run(prog_cmd, CMD = `echo 123 456`)
    run_id_file = run(prog_cmd, CMD = `echo 123 456`, dry_run = true)[2]
    file_info = Pipelines.parse_run_id_file(run_id_file)
    @test isempty(file_info)

    cmd = pipeline(`echo 123 456`, "123456.txt")
    run(prog_cmd, CMD = cmd)
    run_id_file = run(prog_cmd, CMD = cmd, dry_run = true)[2]
    file_info = Pipelines.parse_run_id_file(run_id_file)
    @test haskey(file_info, "123456.txt")

    @test_warn "Skipped finished" run(prog_cmd, CMD = cmd)

    # touching the file, re run
    touch("123456.txt")
    @test_warn "Started: Program Cmd" run(prog_cmd, CMD = cmd)

    cmd = `echo 123456.txt`
    run(prog_cmd, CMD = cmd)
    run_id_file = run(prog_cmd, CMD = cmd, dry_run = true)[2]
    file_info = Pipelines.parse_run_id_file(run_id_file)
    @test haskey(file_info, "123456.txt")

    touch("123.txt")

    cmd = `echo -J=123456.txt,123.txt`
    run(prog_cmd, CMD = cmd)
    run_id_file = run(prog_cmd, CMD = cmd, dry_run = true)[2]
    file_info = Pipelines.parse_run_id_file(run_id_file)
    @test haskey(file_info, "123456.txt")

    cmd = `echo 123456.txt,123.txt`
    run(prog_cmd, CMD = cmd)
    run_id_file = run(prog_cmd, CMD = cmd, dry_run = true)[2]
    file_info = Pipelines.parse_run_id_file(run_id_file)
    @test haskey(file_info, "123456.txt")
    @test haskey(file_info, "123.txt")

    touch("123.txt")
    @test_warn "Started: Program Cmd" run(prog_cmd, CMD = cmd)
    @test_warn "Skipped finished" run(prog_cmd, CMD = cmd)
end