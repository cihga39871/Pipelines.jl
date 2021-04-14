mutable struct JuliaProgram <: Program
	name::String
	id_file::String
	info_before::String
	info_after::String
    cmd_dependencies::Vector{CmdDependency}
	inputs::Vector{String}
	validate_inputs::Function
	prerequisites::Function
    main::Function
	infer_outputs::Function
	outputs::Vector{String}
	validate_outputs::Function
	wrap_up::Function
end

"""
# Struct

	mutable struct JuliaProgram <: Program
		name::String
		id_file::String
		info_before::String
		info_after::String
		cmd_dependencies::Vector{CmdDependency}
		inputs::Vector{String}
		validate_inputs::Function
		prerequisites::Function
		main::Function
		infer_outputs::Function
		outputs::Vector{String}
		validate_outputs::Function
		wrap_up::Function
	end

# Methods

	JuliaProgram(;
		name::String               = "Unnamed Command Program",
		id_file::String            = "",
		info_before::String        = "auto",
		info_after::String         = "auto",
		cmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),
		inputs::Vector{String}     = Vector{String}(),
		validate_inputs::Function  = do_nothing,  # positional arguments: inputs::Dict{String}
		prerequisites::Function    = do_nothing,  # positional arguments: inputs, outputs::Dict{String}
		main::Function             = do_nothing,  # positional arguments: inputs, outputs::Dict{String}
		infer_outputs::Function    = do_nothing,  # positional arguments: inputs::Dict{String}
		outputs::Vector{String}    = Vector{String}(),
		validate_outputs::Function = do_nothing  # positional arguments: outputs::Dict{String},
		wrap_up::Function          = do_nothing  # positional arguments: inputs, outputs::Dict{String}
	) -> JuliaProgram

Julia program template. To run a `JuliaProgram`, use `run(::JuliaProgram; kwargs...).`

# Arguments

- `name::String`: Program name.

- `id_file::String`: The prefix of *run ID file*. To prevent from running the program with the same inputs and outputs twice, it will generate a unique *run ID file* after a successful run.

- `info_before::String`: Print it when starting the program.

- `info_after::String`: Print it when finishing the program.

- `cmd_dependencies::Vector{CmdDependency}`: Any command dependencies used in the program.

- `inputs` and `outputs`: *keys* (`Vector{String}`) of dicts which are required by `run(::JuliaProgram, inputs::Dict{String}, outputs::Dict{String})`. See details below.

  > `JuliaProgram` stores a Julia function, with two arguments `inputs::Dict{String}, outputs::Dict{String}`. The keys of the two arguments should be set in `inputs::Vector{String}` and `outputs::Vector{String}`.

  > Caution: the returned value of `p.main` will be assigned to new `outputs`. Please ensure the returned value is `Dict{String}` with proper keys.

- `validate_inputs::Function`: A function to validate inputs. It takes *one* argument `Dict{String}` whose keys are the same as `inputs`. If validation fail, throw error or return false.

- `prerequisites`: A function to run just before the main command. It prepares necessary things, such as creating directories. It takes *two* arguments `Dict{String}` whose keys are the same as `inputs` and `outputs`, respectively.

- `main::Function`: The main julia function. It takes *two* arguments `Dict{String}` whose keys are the same as `inputs` and `outputs`, respectively.

  > Caution: the returned value of `p.main` will be assigned to new `outputs`. Please ensure the returned value is `Dict{String}` with proper keys.

- `infer_outputs::Function`: A function to infer outputs from inputs. It takes *one* argument `Dict{String}` whose keys are the same as `inputs`.

- `validate_outputs::Function`: A function to validate outputs. It takes *one* argument `Dict{String}` whose keys are the same as `outputs`. If validation fail, throw error or return false.

- `wrap_up::Function`: the last function to run. It takes *two* arguments `Dict{String}` whose keys are the same as `inputs` and `outputs`, respectively.

# Example

	p = JuliaProgram(
		id_file = "id_file",
		inputs = ["a", "b"],
		outputs = ["c"],
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
"""
function JuliaProgram(;
	name::String               = "Julia Program",
	id_file::String            = "",
	info_before::String        = "auto",
	info_after::String         = "auto",
	cmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),
	inputs::Vector{String}     = Vector{String}(),
	validate_inputs::Function  = do_nothing,  # positional arguments: inputs::Dict{String}
	prerequisites::Function    = do_nothing,  # positional arguments: inputs, outputs::Dict{String}
	main::Function             = do_nothing,  # positional arguments: inputs, outputs::Dict{String},
	infer_outputs::Function    = do_nothing,  # positional arguments: inputs::Dict{String}
	outputs::Vector{String}    = Vector{String}(),
	validate_outputs::Function = do_nothing,  # positional arguments: outputs::Dict{String}
	wrap_up::Function          = do_nothing  # positional arguments: inputs, outputs::Dict{String}
)
	JuliaProgram(
		name,
		id_file,
		info_before,
		info_after,
		cmd_dependencies,
		inputs,
		validate_inputs,
		prerequisites,
		main,
		infer_outputs,
		outputs,
		validate_outputs,
		wrap_up
	)
end

"""
	run(
		p::JuliaProgram;
		inputs::Dict{String}=Dict{String, Cmd}(),
		outputs::Dict{String}=Dict{String, Cmd}(),
		skip_when_done::Bool=true,
		check_dependencies::Bool=true,
		stdout=nothing,
		stderr=nothing,
		append::Bool=false,
		verbose::Bool=true,
		touch_run_id_file::Bool=true,
		dry_run::Bool=false
	) -> (success::Bool, outputs::Dict{String})

	run(
		p::JuliaProgram,
		inputs::Dict{String},
		outputs::Dict{String};
		kwargs...
	)

	run(
		p::JuliaProgram,
		inputs::Dict{String},
		kwargs...
	)  # only usable when `p.infer_outputs` is defined.

Run Julia Program (JuliaProgram) using specified `inputs` and `outputs`.

Return `(success::Bool, outputs::Dict{String})`

- `p::JuliaProgram`: the command program template.

- `inputs::Dict{String}` and `outputs::Dict{String}`: `JuliaProgram` stores a Julia function, with two arguments `inputs::Dict{String}, outputs::Dict{String}`. The keys of the two arguments should be the same as `p.inputs::Vector{String}` and `p.outputs::Vector{String}`.

  > Caution: the returned value of `p.main` will be assigned to new `outputs`. Please ensure the returned value of `p.main` is `Dict{String}` with proper keys.

- `skip_when_done::Bool = true`: Skip running the program and return `true` if it has been done before (the `run_id_file` exists and `p.validate_outputs(outputs)` passes.)

- `check_dependencies::Bool = true`: check dependencies for `p` (`p.cmd_dependencies`).

- `stdout`, `stderr` and `append`: Redirect the program output to files, the same behavior as `pipeline(cmd; stdout=stdout, stderr=stderr, append=append)`. Caution: use after checking whether the original command has redirection.

- `verbose::Bool = true`: If `true`, print info and error messages. If `false`, print error messages only.

- `touch_run_id_file::Bool = true`: If `true`, touch a unique run ID file, which indicate the program is successfully run with given inputs and outputs. If `false`, the next time running the program, `skip_when_done=true` will not take effect.

- `dry_run::Bool = false`: do not run the command, return `(fake_outputs::Dict{String}, run_id_file::String)`.

### Workflow

1. Validate compatibility between `p` and `inputs/outputs`.

2. Check whether the program has run before. (`skip_when_done`, `p.validate_outputs(outputs)`)

3. Check command dependencies. (`check_dependencies`, `p.cmd_dependencies`)

4. Validate `inputs`. (`p.validate_inputs(inputs)`)

5. Preparing before running main command. (`p.prerequisites(inputs, outputs)`)

6. Run main function and the returned value will be assigned to new `outputs`. (`outputs = p.main(inputs, outputs)`)

7. Validate `outputs`. (`p.validate_outputs(outputs)`)

8. Wrap up. (`p.wrap_up(inputs, outputs)`)

9. Success, touch run id file, and return true. (`touch_run_id_file::Bool`)

# Example

	p = JuliaProgram(
		id_file = "id_file",
		inputs = ["a", "b"],
		outputs = ["out"],
		main = (inputs, outputs) -> begin
			a = inputs["a"]
			b = inputs["b"]
			println("inputs are ", a, " and ", b)
			println("You can also use info in outputs: ", outputs["out"])

			println("The returned value will be assigned to a new outputs")
			return Dict{String,Any}("out" => b^2)
		end
	)

	inputs = Dict(
		"a" => `in1`,
		"b" => 2
	)

	outputs = Dict(
		"out" => "will_be_replaced"
	)

	success, outputs = run(p, inputs, outputs;
		touch_run_id_file = false
	) # outputs will be refreshed
"""
function Base.run(
	p::JuliaProgram;
	inputs::Dict{String}=Dict{String, Cmd}(),
	outputs::Dict{String}=Dict{String, Cmd}(),
	skip_when_done::Bool=true,
	check_dependencies::Bool=true,
	stdout=nothing,
	stderr=nothing,
	append::Bool=false,
	verbose::Bool=true,
	touch_run_id_file::Bool=true,
	dry_run::Bool=false
)
	# check keyword consistency (keys in inputs and outputs compatible with the function)
	check_keywords(p, inputs, outputs)
	# run id based on inputs and outputs
	run_id = generate_run_uuid(inputs, outputs)
	run_id_file = p.id_file * "." * string(run_id)

	# whether dry run
	if dry_run
		@goto dry_run_start
	end

	if verbose
		if p.info_before == "auto" || p.info_before == ""
			@info timestamp() * "Starting program: $(p.name)" run_id inputs outputs
		else
			@info timestamp() * p.info_before run_id inputs outputs
		end
	end

	# skip when done
	if skip_when_done && isfile(run_id_file) && isok(p.validate_outputs(outputs))
		@warn timestamp() * "Skipping finished program: $(p.name)" run_id inputs outputs
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
		@error timestamp() * "ProgramInputValidationError: $(p.name): fail to validate inputs (before running the main command). See error messages above." validation_function=p.validate_inputs run_id inputs outputs
		error("ProgramInputValidationError")
		return false, outputs
	end

	@label dry_run_start

	if dry_run
		return (outputs, run_id_file)
		# dry_run stops
	end

	# preparation: run function prerequisites before the main command
	try
		isok(p.prerequisites(inputs, outputs)) || error("ProgramPrerequisitesError: $(p.name): the prerequisites function returns false.")
	catch e
		rethrow(e)
		@error timestamp() * "ProgramPrerequisitesError: $(p.name): fail to run the prerequisites function (before running the main command). See error messages above." prerequisites=p.prerequisites run_id inputs outputs
		error("ProgramPrerequisitesError")
		return false, outputs
	end

	# run the main command
	outputs = try
		p.main(inputs, outputs)
	catch e
		rethrow(e)
		@error timestamp() * "ProgramRunningError: $(p.name): fail to run the main command. See error messages above." run_id inputs outputs
		error("ProgramRunningError")
		return false, outputs
	end

	# check type of outputs
	if !isa(outputs, Dict{String})
		@error timestamp() * "ProgramMainReturnValueError: $(p.name): the returned value is not a Dict{String}" run_id inputs outputs
		error("ProgramRunningError")
		return false, outputs
	end
	# check keys in outputs::Dict{String}
	check_outputs_keywords(p, outputs)

	# confirmation: validate outputs
	try
		isok(p.validate_outputs(outputs)) || error("ProgramOutputValidationError: $(p.name): the validation function returns false.")
	catch e
		rethrow(e)
		@error timestamp() * "ProgramOutputValidationError: $(p.name): fail to validate outputs (after running the main command). See error messages above." validation_function=p.validate_outputs run_id inputs outputs
		error("ProgramOutputValidationError")
		return false, outputs
	end

	try
		isok(p.wrap_up(inputs, outputs)) || error("ProgramWrapUpError: $(p.name): the wrap_up function returns false.")
	catch e
		rethrow(e)
		@error timestamp() * "ProgramWrapUpError: $(p.name): fail to run the wrap_up function. See error messages above." wrap_up=p.wrap_up run_id inputs outputs
		error("ProgramWrapUpError")
		return false, outputs
	end

	# touch the run_id_file
	touch_run_id_file && touch(run_id_file)

	if verbose
		if p.info_after == "auto" || p.info_after == ""
			@info timestamp() * "Finishing program: $(p.name)" run_id inputs outputs
		else
			@info timestamp() * p.info_after run_id inputs outputs
		end
	end
	return true, outputs
end
