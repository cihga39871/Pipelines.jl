
using Test
using ScopedStreams

using Pipelines

@testset begin

    include("test_Arg.jl")

    ## keyword interpolation
    allputs = Dict(
        "A" => 123,
        "B" => "<A>",
        "C" => "<B><A>"
    )
    @test Pipelines.keyword_interpolation(allputs) == Dict(
        "A" => 123,
        "B" => "123",
        "C" => "123123"
    )

    allputs = Dict(
        "A" => 123,
        "B" => "<C>",
        "C" => "<B><A>"
    )
    @test_throws ErrorException Pipelines.keyword_interpolation(allputs)

    @test_throws ErrorException Pipelines.check_function_methods(+, (Int, String))

    @test Pipelines.parse_verbose("minimum") === :min
    @test Pipelines.parse_verbose(:max) === :all
    @test Pipelines.parse_verbose(:maximum) === :all
    @test Pipelines.parse_verbose(:full) === :all
    @test Pipelines.parse_verbose(:yes) === :all
    @test Pipelines.parse_verbose(:true) === :all
    @test Pipelines.parse_verbose(:no) === :none
    @test Pipelines.parse_verbose(:nothing) === :none
    @test Pipelines.parse_verbose(:null) === :none
    @test Pipelines.parse_verbose(:false) === :none
    @test Pipelines.parse_verbose(nothing) === :none
    @test @test_warn "Cannot determine" Pipelines.parse_verbose(:unkown) === :all

    @test do_nothing() === nothing
    @test Pipelines.isok(nothing) === true
    @test Pipelines.isok("ok") === true
    @test Pipelines.isok("TRUE") === true
    @test Pipelines.isok("True") === true
    @test Pipelines.isok("") === false
    @test Pipelines.isok("n") === false
    @test Pipelines.isok("false") === false
    @test Pipelines.isok("False") === false
    @test Pipelines.isok("F") === false
    @test Pipelines.isok("FALSE") === false
    @test Pipelines.isok(1234) === true

    @test Pipelines.to_xxput_dict(["A" => 123]) == Pipelines.to_xxput_dict(Dict("A" => 123))
    @test_throws ErrorException Pipelines.to_xxput_dict(33)

    @test Pipelines.str(nothing) == ""
    @test Pipelines.str([23,4]) == "23_4"
    @test Pipelines.str(r"abc") == "abc"

    @test Pipelines.to_cmd(r"abc") == `abc`
    @test Pipelines.to_cmd([r"abc"]) == `abc`
    @test Pipelines.to_cmd(["abc"]) == `abc`

    @test split(`abc def`) == ["abc", "def"]

    @test replaceext("abc.def", ".xyz") == "abc.xyz"
    @test replaceext("abc", ".xyz") == "abc.xyz"
    @test replaceext(r"abc.def.ghi", ".xyz") == "abc.def.xyz"

    @test removeext("abc.def") == "abc"
    @test removeext("abc") == "abc"
    @test removeext(r"abc") == "abc"

    ### cmd dependency

    never_dep = CmdDependency(
        exec = `CmdDependencyExpectedNotToExist`,
        test_args = `--version`,
        validate_success = true,
        validate_stdout = x -> occursin(r"^ABCDEFG", x),
        validate_stderr = do_nothing,
        exit_when_fail = true
    )
    @test_throws ErrorException check_dependency(never_dep; exit_when_fail=true)
    @test !check_dependency(never_dep; exit_when_fail=false)

    julia = CmdDependency(
        exec = Base.julia_cmd(),
        test_args = `--version`,
        validate_success = true,
        validate_stdout = x -> occursin(r"^julia version", x),
        validate_stderr = do_nothing,
        exit_when_fail = true
    )
    @test_nowarn display(julia)
    @test_nowarn Base.show(julia)
    @test_nowarn Base.show(deref(stdout), julia)
    @test_nowarn Base.show(deref(stdout), MIME("text/plain"), julia)

    check_dependency(julia)

    @test `$julia` == Base.julia_cmd()

    ### cmd program    
    
    p = CmdProgram(
        cmd_dependencies = [julia],
        id_file = "id_file",
        inputs = [
            "input",
            "input2" => Int,
            "optional_arg" => 5,
            "optional_arg2" => 0.5 => Number
        ],
        outputs = [
            "output" => "<input>.output"
        ],
        cmd = `echo input input2 optional_arg optional_arg2 output`
    )
    @test_nowarn display(p)

    @test check_dependency(p)
    check_dependency(@__MODULE__)
    check_dependency(@__MODULE__; verbose=false, exit_when_fail=false)
    status_dependency(@__MODULE__; verbose=false, exit_when_fail=false)

    never_dep=nothing
    check_dependency(@__MODULE__)
    check_dependency(@__MODULE__; verbose=false, exit_when_fail=false)
    status_dependency(@__MODULE__; verbose=false, exit_when_fail=false)
    
    @test p.inputs == ["input", "input2", "optional_arg", "optional_arg2"]
    @test p.default_inputs == [nothing, nothing, 5, 0.5]
    @test p.outputs == ["output"]
    @test p.default_outputs == ["<input>.output"]
    @test p.output_types == [Any]

    inputs = Dict(
        "input" => `in1`,
        "input2" => 2
    )
    @test Pipelines.infer_outputs(p, inputs) == Dict(
        "output" => "in1.output"
    )
    @test Pipelines.infer_outputs(p, ["input"=>"in1", "input2"=>4]) == Dict(
        "output" => "in1.output"
    )

    outputs = Dict(
        "output" => "out"
    )

    @test Pipelines.prepare_cmd(p, inputs, outputs) == `echo in1 2 5 0.5 out`
    @test Pipelines.prepare_cmd(p, inputs, "output" => "out") == `echo in1 2 5 0.5 out`

    @test run(p,
        inputs = inputs,
        outputs = outputs,
        skip_when_done = false,
        verbose = true,
        touch_run_id_file = false
    )[1]

    @test prog_run(p,
        inputs,
        outputs;
        skip_when_done = false,
        verbose = false,
        touch_run_id_file = false
    )[1]

    p_nooutput = CmdProgram(
        cmd_dependencies = [julia],
        id_file = "id_file",
        inputs = ["input", "input2", "optional_arg", "optional_arg2"],
        cmd = `echo input input2`
    )

    @test_throws ErrorException run(p_nooutput,
        inputs;
        skip_when_done = false,
        verbose = false,
        touch_run_id_file = false
    )

    cmd, run_id_file = run(p,
        inputs = Dict(
            "input" => `in1`,
            "input2" => 5.0
        ),
        outputs = Dict(
            "output" => `out`
        ),
        dry_run=true
    )
    rm(run_id_file, force=true)

    program_bowtie2 = CmdProgram(
        name = "Bowtie2 Mapping",
        id_file = ".bowtie2",

        inputs = ["FASTQ", "REF"],
        validate_inputs = inputs -> begin
            check_dependency_file(inputs["FASTQ"]) &&
            check_dependency_file(inputs["REF"])
        end,

        prerequisites = (inputs, outputs) -> begin
            mkpath(dirname(to_str(outputs["BAM"])))
        end,

        outputs = ["BAM"],
        infer_outputs = inputs -> begin
            Dict("BAM" => str(inputs["FASTQ"]) * ".bam")
        end,
        validate_outputs = outputs -> begin
            check_dependency_file(outputs["BAM"])
        end,

        cmd = pipeline(`bowtie2 -x REF -q FASTQ`, `samtools sort -O bam -o BAM`),

        wrap_up = (inputs, outputs) -> run(`samtools index $(outputs["BAM"])`)
    )


    ### julia program
    p = JuliaProgram(
        id_file = "id_file",
        inputs = [
            "a",
            "b" => Int
        ],
        outputs = [
            "c" => "<a>.<b>"
        ],
        main = (inputs, outputs) -> begin
            a = inputs["a"]
            b = inputs["b"]
            println("inputs are ", a, " and ", b)
            println("You can also use info in outputs: ", outputs["c"])
            println("The returned value will be assigned to a new outputs")

            return Dict{String,Any}("c" => b^2)
        end
    )

    @test_nowarn display(p)
    @test_nowarn Base.show(p)
    @test_nowarn Base.show(deref(stdout), p)
    @test_nowarn Base.show(deref(stdout), MIME("text/plain"), p)


    inputs = Dict(
        "a" => `in1`,
        "b" => 2
    )

    outputs = Dict(
        "c" => "out"
    )

    success, outputs = run(p, inputs, outputs;
        touch_run_id_file = false
    ) # outputs will be refreshed
    @test success

    success, outputs = run(p, inputs;
        touch_run_id_file = false
    ) # outputs will be refreshed
    @test success


    ### run program without inputs outputs
    p = JuliaProgram(
        id_file = "id_file",
        inputs = [
            "a" => 10.6,
            "b" =>  5 => Int
        ],
        outputs = "c" => "<a>.<b>",
        main = (inputs, outputs) -> begin
            return Dict{String,Any}("c" => inputs["b"]^2)
        end
    )
    success, outputs = run(p;
        touch_run_id_file = false
    )
    @test success

    inputs, outputs = Pipelines.xxputs_completion_and_check(p, Dict{String, Any}(), Dict{String, Any}())

    p = JuliaProgram(
        id_file = "id_file",
        inputs = [
            "a" => 10.6 => Float64,
            "b" =>  5 => Int
        ],
        outputs = "c" => "<a>.<b>",
        main = (inputs, outputs) -> begin
            return Dict{String,Any}("c" => inputs["b"]^2)
        end
    )
    success, outputs = run(p;
        touch_run_id_file = false
    )
    @test success

    success, outputs = run(p,
        inputs = Dict(
            "a" => 10.6,
            "b" =>  2.0
        ),
        touch_run_id_file = false
    )
    @test success

    success, outputs = run(p,
        inputs = :a => 8,
        touch_run_id_file = false
    )
    @test success

    ### testing dir, stdout/err/log

    ## julia program
    p = JuliaProgram(
        id_file = "id_file",
        inputs = [
            :a => 10.6 => Float64,
            :b =>  5 => Int
        ],
        outputs = "c" => "<a>.<b>",
        main = (inputs, outputs) -> begin
            println(stdout, "stdout> ", inputs["a"])
            println(stderr, "stderr> ", inputs["b"])
            @info outputs
            return Dict{String,Any}("c" => inputs["b"]^2)
        end
    )

    tmp = mktempdir()
    cd(tmp)
    working_dir = mktempdir()
    @test run(p, touch_run_id_file=false, stdout="out.txt", stdlog="log.txt")[1]

    @test read("out.txt", String) == "stdout> 10.6\n"
    @test length(read("log.txt", String)) > 500

    @test run(p, "a" => 6.6, touch_run_id_file=false, stderr = "err.txt", dir=working_dir)[1]

    @test isfile("err.txt")
    rm("err.txt")

    Pipelines.auto_change_directory(true)
    @test run(p, "a" => 6.6, touch_run_id_file=false, stderr = "err.txt", dir=working_dir)[1]

    # redirect files at working dir
    @test !isfile("err.txt")
    @test isfile(joinpath(working_dir, "err.txt"))
    @test length(read(joinpath(working_dir, "err.txt"), String)) > 500

    Pipelines.auto_change_directory(false)

    # skip when done
    run(p; skip_when_done=true)[1]
    @test_warn "Skipped finished program" run(p; skip_when_done=true)[1]
    @test_warn "Skipped finished program" run(p; skip_when_done=true, verbose=:min)[1]


    ## cmd program
    p = CmdProgram(
        id_file = "id_file",
        inputs = [
            :a => 10.6 => Float64,
            :b =>  5 => Int
        ],
        outputs = "c" => "<a>.<b>",
        cmd = `echo this is stdout a b` & pipeline(`julia -e '@info "this is stderr"'`)
    )
    out_io = open("out.txt", "w+")

    @test run(p, touch_run_id_file=false, stdout=out_io, stderr="err.txt", stdlog="log.txt")[1]

    close(out_io)

    @test read("out.txt", String) == "this is stdout 10.6 5\n"
    @test read("err.txt", String) == "[ Info: this is stderr\n"
    @test !occursin("[ Info: this is stderr", read("log.txt", String))

    rm("out.txt", force=true)
    rm("err.txt", force=true)
    rm("log.txt", force=true)

    p2 = CmdProgram(
        id_file = "id_file",
        inputs = [
            :a => 10.6 => Float64,
            :b =>  5 => Int
        ],
        outputs = "c" => "<a>.<b>",
        cmd = pipeline(`echo 123`, "out2.txt")
    )
    @test run(p2, touch_run_id_file=false, stdout=stdout)[1]  # stdout=stdout will not take effect because redirection was done in out2.txt
    @test isfile("out2.txt") && read("out2.txt", String) == "123\n"
    rm("out2.txt", force=true)


    ## version 0.6

    # test retry
    p_error = CmdProgram(
        id_file = "id_file",
        cmd = `julia --abcdefg`
    )
    @test_throws ProcessFailedException run(p_error, retry=1, verbose=:min)

    pj_error = JuliaProgram(
        id_file = "id_file",
        main = (x, y) -> begin
            if !isfile("x")
                touch("x")
                error("x not exist, so created")
            end
            return y
        end
    )
    rm("x", force=true)
    res = @test_warn "failed but will retry." run(pj_error, retry=1, verbose=:min, touch_run_id_file=false)
    @test isfile("x")
    rm("x", force=true)


    include("test_v0.8.jl")

    include("test_v0.9.jl")

    # clean up
    cd(homedir())
    rm(tmp, recursive=true)
    rm(working_dir, recursive=true)
    @test stdout isa ScopedStreams.ScopedStream
    @test stderr isa ScopedStreams.ScopedStream
end
