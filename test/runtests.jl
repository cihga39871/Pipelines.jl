
include("../src/Pipelines.jl")

using .Pipelines
using Test

### util
@test Pipelines.parse_default_element("x") == ("x", Any, nothing)
@test Pipelines.parse_default_element("x" => Int) == ("x", Int, nothing)
@test Pipelines.parse_default_element("x" => 5) == ("x", Any, 5)
@test Pipelines.parse_default_element("x" => 5 => Int) == ("x", Int, 5)
@test Pipelines.parse_default_element("x" => Int => 5) == ("x", Int, 5)
@test Pipelines.parse_default_element("x" => Int => DataType) == ("x", DataType, Int)
@test Pipelines.parse_default_element("x" => DataType => Int) == ("x", DataType, Int)
@test Pipelines.parse_default_element( :b => (5.0 => Int)) == ("b", Int, 5)

@test_throws ErrorException Pipelines.parse_default_element("x" => Int => "5")
@test_throws ErrorException Pipelines.parse_default_element("x" => "5" => Int)
@test_throws ErrorException Pipelines.parse_default_element(5)

v = ["A", "B" => Int, "C" => 9.0 => Float64]
@test Pipelines.parse_default(v) == (["A", "B", "C"], DataType[Any, Int64, Float64], Any[nothing, nothing, 9.0])

v = ["A", :B => Int, "C" => 9.0 => Float64]
@test Pipelines.parse_default(v) == (["A", "B", "C"], DataType[Any, Int64, Float64], Any[nothing, nothing, 9.0])

v = ["A", :B => Int => 5.0, "C" => 9.0 => Float64]
@test Pipelines.parse_default(v) == (["A", "B", "C"], DataType[Any, Int64, Float64], Any[nothing, 5, 9.0])

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

### cmd dependency

julia = CmdDependency(
    exec = Base.julia_cmd(),
    test_args = `--version`,
    validate_success = true,
    validate_stdout = x -> occursin(r"^julia version", x),
    validate_stderr = do_nothing,
    exit_when_fail = true
)

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

inputs = Dict(
    "input" => `in1`,
    "input2" => 2
)

outputs = Dict(
    "output" => "out"
)

run(p,
    inputs = inputs,
    outputs = outputs,
    skip_when_done = false,
    verbose = true,
    touch_run_id_file = false
)

run(p,
    inputs,
    outputs;
    skip_when_done = false,
    verbose = false,
    touch_run_id_file = false
)

p_nooutput = CmdProgram(
    cmd_dependencies = [julia],
    id_file = "id_file",
    inputs = ["input", "input2", "optional_arg", "optional_arg2"],
    cmd = `echo input input2`
)

run(p_nooutput,
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

success, outputs = run(p, inputs;
	touch_run_id_file = false
) # outputs will be refreshed


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
success, outputs = run(p,
	inputs = Dict(
		"a" => 10.6,
		"b" =>  2.0
	),
	touch_run_id_file = false
)

success, outputs = run(p,
	inputs = :a => 8,
	touch_run_id_file = false
)

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
run(p, touch_run_id_file=false, stdout="out.txt", stdlog="log.txt")

@test read("out.txt", String) == "stdout> 10.6\n"
@test length(read("log.txt", String)) > 600

run(p, "a" => 6.6, touch_run_id_file=false, stderr = "err.txt", dir=working_dir)

# redirect files at working dir
@test !isfile("err.txt")
@test isfile(joinpath(working_dir, "err.txt"))
@test length(read(joinpath(working_dir, "err.txt"), String)) > 600

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

run(p, touch_run_id_file=false, stdout=out_io, stderr="err.txt", stdlog="log.txt")

close(out_io)

@test read("out.txt", String) == "this is stdout 10.6 5\n"
@test read("err.txt", String) == "[ Info: this is stderr\n"
@test !occursin("[ Info: this is stderr", read("log.txt", String))


## version 0.6

# test retry
p_error = CmdProgram(
	id_file = "id_file",
	cmd = `julia --abcdefg`
)
p_error_res = run(p_error, retry=1, verbose=:min)
@test p_error_res isa Pipelines.StackTraceVector

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
pj_error_res = run(pj_error, retry=1, verbose=:min)
rm("x", force=true)
@test pj_error_res[1]


# clean up
cd(homedir())
rm(tmp, recursive=true)
rm(working_dir, recursive=true)
