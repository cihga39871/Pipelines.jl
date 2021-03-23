
include("../src/Pipelines.jl")

using .Pipelines

julia = CmdDependency(
    exec = Base.julia_cmd(),
    test_args = `--version`,
    validate_success = true,
    validate_stdout = x -> occursin(r"^julia version", x),
    validate_stderr = do_nothing,
    exit_when_fail = true
)

check_dependency(julia)

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
    verbose = false,
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
