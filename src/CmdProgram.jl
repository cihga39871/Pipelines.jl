mutable struct CmdProgram <: Program
	name::String
	id_file::String
	info_before::String
	info_after::String
    cmd_dependencies::Vector{CmdDependency}
	inputs::Vector{String}
	input_types::Vector{DataType}
	default_inputs::Vector
	validate_inputs::Function
	prerequisites::Function
    cmd::Base.AbstractCmd
	infer_outputs::Function
	outputs::Vector{String}
	output_types::Vector{DataType}
	default_outputs::Vector
	validate_outputs::Function
	wrap_up::Function
end

"""
	CmdProgram <: Program

	CmdProgram(;
		name::String               = "Command Program",
		id_file::String            = "",
		info_before::String        = "auto",
		info_after::String         = "auto",
		cmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),
		inputs                     = Vector{String}(),
		validate_inputs::Function  = do_nothing,  # positional arguments: inputs::Dict{String, ValidInputTypes}
		prerequisites::Function    = do_nothing,  # positional arguments: inputs, outputs::Dict{String, ValidInputTypes}
		cmd::Base.AbstractCmd      = ``,
		infer_outputs::Function    = do_nothing,  # positional arguments: inputs::Dict{String, ValidInputTypes}
		outputs                    = Vector{String}(),
		validate_outputs::Function = do_nothing  # positional arguments: outputs::Dict{String, ValidInputTypes},
		wrap_up::Function          = do_nothing  # positional arguments: inputs, outputs::Dict{String, ValidInputTypes}
	) -> CmdProgram

Command program template. To run a `CmdProgram`, use `run(::CmdProgram; kwargs...).`

# Arguments

- `name::String`: Program name.

- `id_file::String`: The prefix of *run ID file*. To prevent from running the program with the same inputs and outputs twice, it will generate a unique *run ID file* after a successful run.

- `info_before::String`: Print it when starting the program.

- `info_after::String`: Print it when finishing the program.

- `cmd_dependencies::Vector{CmdDependency}`: Any command dependencies used in the program.

- `inputs` and `outputs`: Elements (or vectors containing elements) in the following format: (1) `keyword` (2) `keyword => data_type` (3) `keyword => default_value` (4) `keyword => default_value => data_type`.

  `keyword` is an argument name, can be `String` or `Symbol`.

  `default_value` is optional. If set, users may not provide this argument when running. Elsewise, users have to provide it. Caution: `nothing` is preserved and means default value not set. If `String`, it can contain other keywords, but need to quote using '<>', such as `"<arg>.txt"`

  `data_type` is optional. If set, the value provided have to be this data type, or an error will throw.

  *HOW DOES THIS WORK?*

  > `CmdProgram` stores a command template. In the template, replaceable portions are occupied by *keywords*, and all keywords are set in `inputs` and `outputs`.
  > `keyword`s will be replaced before running the program. Users need to provide a dictionary of `keyword::String => value` in `run(::Program, inputs::Dict{String}, outputs::Dict{String})`.

- `validate_inputs::Function`: A function to validate inputs. It takes *one* argument `Dict{String, ValidInputTypes}` whose keys are the same as `inputs`. If validation fail, throw error or return false.

- `prerequisites`: A function to run just before the main command. It prepares necessary things, such as creating directories. It takes *two* arguments `Dict{String, ValidInputTypes}` whose keys are the same as `inputs` and `outputs`, respectively.

- `cmd::AbstractCmd`: The main command template. In the template, keywords in `inputs::Vector{String}` and `outputs::Vector{String}` will be replaced when envoking `run(::CmdProgram, inputs::Dict{String, ValidInputTypes}, outputs::Dict{String, ValidInputTypes})`.

- `infer_outputs::Function`: A function to infer outputs from inputs. It takes *one* argument `Dict{String, ValidInputTypes}` whose keys are the same as `inputs`.

- `validate_outputs::Function`: A function to validate outputs. It takes *one* argument `Dict{String, ValidInputTypes}` whose keys are the same as `outputs`. If validation fail, throw error or return false.

- `wrap_up::Function`: the last function to run. It takes *two* arguments `Dict{String, ValidInputTypes}` whose keys are the same as `inputs` and `outputs`, respectively.

# Example

	p = CmdProgram(
		id_file = "id_file",
		inputs = [
			"input",
			"input2" => Int,
			"optional_arg" => 5,
			"optional_arg2" => 0.5 => Number
		],
		outputs =
			"output" => "<input>.output"
		,
		cmd = `echo input input2 optional_arg optional_arg2 output`
	)

	inputs = Dict(
		"input" => `in1`,
		"input2" => 2
	)

	outputs = Dict(
		"output" => "out"
	)

	run(p, inputs, outputs;
		touch_run_id_file = false
	)
"""
function CmdProgram(;
	name::String               = "Command Program",
	id_file::String            = "",
	info_before::String        = "auto",
	info_after::String         = "auto",
	cmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),
	inputs                     = Vector{String}(),
	validate_inputs::Function  = do_nothing,  # positional arguments: inputs::Dict{String}
	prerequisites::Function    = do_nothing,  # positional arguments: inputs, outputs::Dict{String}
	cmd::Base.AbstractCmd      = ``,
	infer_outputs::Function    = do_nothing,  # positional arguments: inputs::Dict{String}
	outputs                    = Vector{String}(),
	validate_outputs::Function = do_nothing,  # positional arguments: outputs::Dict{String}
	wrap_up::Function          = do_nothing  # positional arguments: inputs, outputs::Dict{String}
)
	inputs, input_types, default_inputs = parse_default(inputs)
	outputs, output_types, default_outputs = parse_default(outputs)

	CmdProgram(
		name,
		id_file,
		info_before,
		info_after,
		cmd_dependencies,
		inputs,
		input_types,
		default_inputs,
		validate_inputs,
		prerequisites,
		cmd,
		infer_outputs,
		outputs,
		output_types,
		default_outputs,
		validate_outputs,
		wrap_up
	)
end

"""
	run(
		p::CmdProgram;
		inputs=Dict{String}(),
		outputs=Dict{String}(),
		skip_when_done::Bool=true,
		check_dependencies::Bool=true,
		stdout=nothing,
		stderr=nothing,
		append::Bool=false,
		verbose::Bool=true,
		touch_run_id_file::Bool=true,
		dry_run::Bool=false
	) -> (success::Bool, outputs::Dict{String})

	run(p::CmdProgram, inputs, outputs; kwargs...)

	run(p::CmdProgram, inputs; kwargs...)
	)  # only usable when `p.infer_outputs` is defined, or default outputs are set in `p`.

Run Command Program (CmdProgram) using specified `inputs` and `outputs`.

Return `(success::Bool, outputs::Dict{String})`

- `p::CmdProgram`: the command program template.

- `inputs::Dict{String}` and `outputs::Dict{String}`: `p::CmdProgram` stores a command template. In the template, replaceable portions are occupied by *keywords*, and all keywords can be found at `p.inputs` and `p.outputs` string vectors. Here, `inputs` and `outputs` are `Dict(keyword::String => replacement)`. The replacements do not have a length limit, unless a *keyword* refers to a filename (length == 1).

  > If data types of `inputs` and `outputs` are not `Dict{String}`, they will be converted as far as possible. If the conversion fails, program will throw an error.

- `skip_when_done::Bool = true`: Skip running the program and return `true` if it has been done before (the `run_id_file` exists and `p.validate_outputs(outputs)` passes.)

- `check_dependencies::Bool = true`: check dependencies for `p` (`p.cmd_dependencies`).

- `stdout`, `stderr` and `append`: Redirect the program output to files, the same behavior as `pipeline(cmd; stdout=stdout, stderr=stderr, append=append)`. Caution: use after checking whether the original command has redirection.

- `verbose::Bool = true`: If `true`, print info and error messages. If `false`, print error messages only.

- `touch_run_id_file::Bool = true`: If `true`, touch a unique run ID file, which indicate the program is successfully run with given inputs and outputs. If `false`, the next time running the program, `skip_when_done=true` will not take effect.

- `dry_run::Bool = false`: do not run the command, return `(command::AbstractCmd, run_id_file::String)`.

### Workflow

1. Validate compatibility between `p` and `inputs/outputs`.

2. Check whether the program has run before. (`skip_when_done`, `p.validate_outputs(outputs)`)

3. Check command dependencies. (`check_dependencies`, `p.cmd_dependencies`)

4. Validate `inputs`. (`p.validate_inputs(inputs)`)

5. Generate runnable command from `p` and `inputs/outputs`. (`stdout`, `stderr`, `append`)

6. Preparing before running main command. (`p.prerequisites(inputs, outputs)`)

7. Run command generated in #5.

8. Validate `outputs`. (`p.validate_outputs(outputs)`)

9. Wrap up. (`p.wrap_up(inputs, outputs)`)

10. Success, touch run id file, and return true. (`touch_run_id_file::Bool`)

# Example

	p = CmdProgram(
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

	run(p, inputs, outputs;
		touch_run_id_file = false
	)
"""
function Base.run(
	p::CmdProgram;
	inputs=Dict{String, Any}(),
	outputs=Dict{String, Any}(),
	skip_when_done::Bool=true,
	check_dependencies::Bool=true,
	stdout=nothing,
	stderr=nothing,
	append::Bool=false,
	verbose::Bool=true,
	touch_run_id_file::Bool=true,
	dry_run::Bool=false
)
	# input/output completion
	inputs, outputs = xxputs_completion_and_check(p, inputs, outputs)

	# run id based on inputs and outputs
	run_id = generate_run_uuid(inputs, outputs)
	run_id_file = p.id_file * "." * string(run_id)

	# whether dry run
	if dry_run
		@goto dry_run_start
	end

	if verbose
		if p.info_before == "auto" || p.info_before == ""
			@info timestamp() * "Starting program: $(p.name)" command_template=p.cmd run_id inputs outputs
		else
			@info timestamp() * p.info_before command_template=p.cmd run_id inputs outputs
		end
	end

	# skip when done
	if skip_when_done && isfile(run_id_file) && isok(p.validate_outputs(outputs))
		@warn timestamp() * "Skipping finished program: $(p.name)" command_template=p.cmd run_id inputs outputs
		return true, outputs
	end

	# check dependencies
	if check_dependencies
		foreach(check_dependency, p.cmd_dependencies)
	end

	# preparation: remove run id file
	isfile(run_id_file) && rm(run_id_file)

	# preparation: validate inputs
	try
		isok(p.validate_inputs(inputs)) || error("ProgramInputValidationError: $(p.name): the prerequisites function returns false.")
	catch e
		rethrow(e)
		@error timestamp() * "ProgramInputValidationError: $(p.name): fail to validate inputs (before running the main command). See error messages above." validation_function=p.validate_inputs command_template=p.cmd run_id inputs outputs
		error("ProgramInputValidationError")
		return false, outputs
	end

	@label dry_run_start
	# preparation: replace inputs and outputs in cmd, including redirecting files
	cmd = prepare_cmd(p, inputs, outputs)
	# redirecting to stdout/stderr if specified.
	if !(isnothing(stdout) && isnothing(stderr))
		cmd = pipeline(cmd, stdout=stdout, stderr=stderr, append=append)
	end

	if dry_run
		return (cmd, run_id_file)
		# dry_run stops
	end

	# preparation: run function prerequisites before the main command
	try
		isok(p.prerequisites(inputs, outputs)) || error("ProgramPrerequisitesError: $(p.name): the prerequisites function returns false.")
	catch e
		rethrow(e)
		@error timestamp() * "ProgramPrerequisitesError: $(p.name): fail to run the prerequisites function (before running the main command). See error messages above." prerequisites=p.prerequisites command_template=p.cmd command_running=cmd run_id inputs outputs
		error("ProgramPrerequisitesError")
		return false, outputs
	end

	# run the main command
	try
		run(cmd)
	catch e
		rethrow(e)
		@error timestamp() * "ProgramRunningError: $(p.name): fail to run the main command. See error messages above." prerequisites=p.prerequisites command_template=p.cmd command_running=cmd run_id inputs outputs
		error("ProgramRunningError")
		return false, outputs
	end

	# confirmation: validate outputs
	try
		isok(p.validate_outputs(outputs)) || error("ProgramOutputValidationError: $(p.name): the validation function returns false.")
	catch e
		rethrow(e)
		@error timestamp() * "ProgramOutputValidationError: $(p.name): fail to validate outputs (after running the main command). See error messages above." validation_function=p.validate_outputs command_template=p.cmd command_running=cmd run_id inputs outputs
		error("ProgramOutputValidationError")
		return false, outputs
	end

	try
		isok(p.wrap_up(inputs, outputs)) || error("ProgramWrapUpError: $(p.name): the wrap_up function returns false.")
	catch e
		rethrow(e)
		@error timestamp() * "ProgramWrapUpError: $(p.name): fail to run the wrap_up function. See error messages above." wrap_up=p.wrap_up command_template=p.cmd command_running=cmd run_id inputs outputs
		error("ProgramWrapUpError")
		return false, outputs
	end

	# touch the run_id_file
	touch_run_id_file && touch(run_id_file)

	if verbose
		if p.info_after == "auto" || p.info_after == ""
			@info timestamp() * "Finishing program: $(p.name)" command_template=p.cmd command_running=cmd run_id inputs outputs
		else
			@info timestamp() * p.info_after command_template=p.cmd command_running=cmd run_id inputs outputs
		end
	end
	return true, outputs
end

function prepare_cmd(c::Cmd, inputs::Dict{String}, outputs::Dict{String})
	cmd = deepcopy(c)
	n = length(cmd.exec)
	for i = n:-1:1
		replacement = get(inputs, cmd.exec[i], nothing) |> to_cmd
		if isnothing(replacement)
			replacement = get(outputs, cmd.exec[i], nothing) |> to_cmd
			if !isnothing(replacement)
				splice!(cmd.exec, i, replacement)
			end
		else
			splice!(cmd.exec, i, replacement)
		end
	end
	cmd
end

function prepare_cmd(c::Base.CmdRedirect, inputs::Dict{String}, outputs::Dict{String})
	Base.CmdRedirect(
		prepare_cmd(c.cmd, inputs, outputs),
		prepare_cmd(c.handle, inputs, outputs),
		c.stream_no,
		c.readable
	)
end
function prepare_cmd(c::T, inputs::Dict{String}, outputs::Dict{String}) where T <: Union{Base.OrCmds, Base.ErrOrCmds, Base.AndCmds}
	T(
		prepare_cmd(c.a, inputs, outputs),
		prepare_cmd(c.b, inputs, outputs)
	)
end
function prepare_cmd(h::Base.FileRedirect, inputs::Dict{String}, outputs::Dict{String})
	replacement = get(outputs, h.filename, nothing) |> to_cmd
	if isnothing(replacement)
		replacement = get(inputs, h.filename, nothing) |> to_cmd
		if !isnothing(replacement)
			@goto to_replace
		end
	else
		@goto to_replace
	end
	return h  # nothing change

	@label to_replace
	if length(replacement.exec) == 1
		return Base.FileRedirect(replacement.exec[1], h.append)
	else
		error("ProgramRedirectError: a key ($(h.filename)) in inputs or outputs found in redirecting filename, but the provided value contain multiple elements. Please only provide one element. For example, quote the filename when contain special characters.")
	end
end

"""
	prepare_cmd(p::CmdProgram, inputs, outputs)

Prepare the runable command. Keywords in CmdProgram will be given to values of inputs/outputs.
"""
function prepare_cmd(p::CmdProgram, inputs::Dict{String}, outputs::Dict{String})
	prepare_cmd(p.cmd, inputs, outputs)
end

function prepare_cmd(p::CmdProgram, inputs, outputs)
	prepare_cmd(p.cmd, to_xxput_dict(inputs), to_xxput_dict(outputs))
end
