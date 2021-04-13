
include("../src/Pipelines.jl")

using .Pipelines
using Test

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
    inputs = ["input", "input2"],
    outputs = ["output"],
    cmd = `echo input input2 output`
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
    inputs = ["input", "input2"],
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
        "input2" => `in2`
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
	cmd_dependencies = [julia],
	id_file = "id_file",
	inputs = ["a", "b"],
	outputs = ["c"],
	main = (inputs, outputs) -> begin
		println("inputs are ", inputs["a"], " and ", inputs["b"])
		println("output is ", outputs["c"])
	end
)

inputs = Dict(
	"a" => `in1`,
	"b" => 2
)

outputs = Dict(
	"c" => "out"
)

run(p, inputs, outputs;
	touch_run_id_file = false
)
