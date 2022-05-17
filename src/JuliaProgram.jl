mutable struct JuliaProgram <: Program
	name::String
	id_file::String
	info_before::String
	info_after::String
    cmd_dependencies::Vector{CmdDependency}
	inputs::Vector{String}
	input_types::Vector{Type}
	default_inputs::Vector
	validate_inputs::Function
	prerequisites::Function
    main::Function
	infer_outputs::Function
	outputs::Vector{String}
	output_types::Vector{Type}
	default_outputs::Vector
	validate_outputs::Function
	wrap_up::Function

	function JuliaProgram(name, id_file, info_before, info_after, cmd_dependencies, inputs, input_types, default_inputs, validate_inputs, prerequisites, main, infer_outputs, outputs, output_types, default_outputs, validate_outputs, wrap_up)

		check_reserved_xxputs(inputs)
		check_reserved_xxputs(outputs)

		check_function_methods(validate_inputs, (Dict, ), "validate_inputs")
		check_function_methods(prerequisites, (Dict, Dict), "prerequisites")
		check_function_methods(main, (Dict, Dict), "main")
		check_function_methods(infer_outputs, (Dict, ), "infer_outputs")
		check_function_methods(validate_outputs, (Dict, ), "validate_outputs")
		check_function_methods(wrap_up, (Dict, Dict), "wrap_up")

		new(name, id_file, info_before, info_after, cmd_dependencies, inputs, input_types, default_inputs, validate_inputs, prerequisites, main, infer_outputs, outputs, output_types, default_outputs, validate_outputs, wrap_up)
	end
end

"""
	JuliaProgram <: Program

	JuliaProgram(;
		name::String               = "Unnamed Command Program",
		id_file::String            = "",
		info_before::String        = "auto",
		info_after::String         = "auto",
		cmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),
		inputs                     = Vector{String}(),
		validate_inputs::Function  = do_nothing,  # positional arguments: inputs::Dict{String}
		prerequisites::Function    = do_nothing,  # positional arguments: inputs, outputs::Dict{String}
		main::Function             = do_nothing,  # positional arguments: inputs, outputs::Dict{String}
		infer_outputs::Function    = do_nothing,  # positional arguments: inputs::Dict{String}
		outputs                    = Vector{String}(),
		validate_outputs::Function = do_nothing,  # positional arguments: outputs::Dict{String}
		wrap_up::Function          = do_nothing  # positional arguments: inputs, outputs::Dict{String}
	) -> JuliaProgram

Julia program template. To run a `JuliaProgram`, use `run(::JuliaProgram; kwargs...).`

# Arguments

- `name::String`: Program name.

- `id_file::String`: The prefix of *run ID file*. To prevent from running the program with the same inputs and outputs twice, it will generate a unique *run ID file* after a successful run.

- `info_before::String`: Print it when the program is started.

- `info_after::String`: Print it when the program is finished.

- `cmd_dependencies::Vector{CmdDependency}`: Any command dependencies used in the program.

- `inputs` and `outputs`: Elements (or vectors containing elements) in the following format: (1) `keyword` (2) `keyword => data_type` (3) `keyword => default_value` (4) `keyword => default_value => data_type`.

  `keyword` is an argument name, can be `String` or `Symbol`.

  `default_value` is optional. If set, users may not provide this argument when running. Elsewise, users have to provide it. Caution: `nothing` is preserved and means default value not set. If `String`, it can contain other keywords, but need to quote using '<>', such as `"<arg>.txt"`

  `data_type` is optional. If set, the value provided have to be this data type, or an error will throw.

  *HOW DOES THIS WORK?*

  > `JuliaProgram` stores a Julia function, with two arguments `inputs::Dict{String}, outputs::Dict{String}`. The keys of the two arguments should be set in `inputs::Vector{String}` and `outputs::Vector{String}`.
  > Caution: the returned value of `p.main` will be assigned to new `outputs` when the returned value is `Dict{String}` with proper keys. Otherwise, a warning will show and the returned value will be the original `outputs::Dict`.

- `validate_inputs::Function`: A function to validate inputs. It takes *one* argument `Dict{String}` whose keys are the same as `inputs`. If validation fail, throw error or return false.

- `prerequisites`: A function to run just before the main command. It prepares necessary things, such as creating directories. It takes *two* arguments `Dict{String}` whose keys are the same as `inputs` and `outputs`, respectively.

- `main::Function`: The main julia function. It takes *two* arguments `Dict{String}` whose keys are the same as `inputs` and `outputs`, respectively.

  > Caution: the returned value of `p.main` will be assigned to new `outputs`. Please ensure the returned value is `Dict{String}` with proper keys.

- `infer_outputs::Function`: A function to infer outputs from inputs. It takes *one* argument `Dict{String}` whose keys are the same as `inputs`.

- `validate_outputs::Function`: A function to validate outputs. It takes *one* argument `Dict{String}` whose keys are the same as `outputs`. If validation fail, throw error or return false.

  > Caution: the returned value of `p.main` will be assigned to new `outputs` when the returned value is `Dict{String}` with proper keys. Otherwise, a warning will show and the returned value will be the original `outputs::Dict`.

- `wrap_up::Function`: the last function to run. It takes *two* arguments `Dict{String}` whose keys are the same as `inputs` and `outputs`, respectively.

# Example

	p = JuliaProgram(
		id_file = "id_file",
		inputs = ["a",
		          "b" => Int],
		outputs = "c" => "<a>.<b>",
		main = (inputs, outputs) -> begin
			a = inputs["a"]
			b = inputs["b"]
			println("inputs are ", a, " and ", b)
			println("You can also use info in outputs: ", outputs["c"])
	        println("The returned value will be assigned to a new outputs")
	        return Dict{String,Any}("c" => b^2)
		end)

	# running the program using `prog_run`: keyword arguments include keys of inputs and outputs
	success, new_out = prog_run(p; a = `in1`, b = 2, c = "out", touch_run_id_file = false)

	@assert new_out != outputs  # outputs will change to the returned value of main function, if the returned value is a Dict and pass `p.validate_outputs`

	# an old way to `run` program: need to create Dicts of inputs and outputs first.
	inputs = Dict("a" => `in1`, "b" => 2)
	outputs = "c" => "out"
	success, new_out = run(p, inputs, outputs; touch_run_id_file = false)
"""
function JuliaProgram(;
	name::String               = "Julia Program",
	id_file::String            = "",
	info_before::String        = "auto",
	info_after::String         = "auto",
	cmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),
	inputs                     = Vector{String}(),
	validate_inputs::Function  = do_nothing,  # positional arguments: inputs::Dict{String}
	prerequisites::Function    = do_nothing,  # positional arguments: inputs, outputs::Dict{String}
	main::Function             = do_nothing,  # positional arguments: inputs, outputs::Dict{String},
	infer_outputs::Function    = do_nothing,  # positional arguments: inputs::Dict{String}
	outputs                    = Vector{String}(),
	validate_outputs::Function = do_nothing,  # positional arguments: outputs::Dict{String}
	wrap_up::Function          = do_nothing  # positional arguments: inputs, outputs::Dict{String}
)
	inputs, input_types, default_inputs = parse_default(inputs)
	outputs, output_types, default_outputs = parse_default(outputs)

	JuliaProgram(
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
		main,
		infer_outputs,
		outputs,
		output_types,
		default_outputs,
		validate_outputs,
		wrap_up
	)
end

function _run(
	p::JuliaProgram;
	inputs=Dict{String, Any}(),
	outputs=Dict{String, Any}(),
	skip_when_done::Bool=true,
	check_dependencies::Bool=true,
	verbose::Union{Bool, Symbol, AbstractString}=true,
	touch_run_id_file::Bool=true,
	dry_run::Bool=false,
	dir::AbstractString=""
)
	# input/output completion
	inputs, outputs = xxputs_completion_and_check(p, inputs, outputs)

	# run id based on inputs and outputs
	run_id = generate_run_uuid(inputs, outputs)
	run_id_file = joinpath(dir, p.id_file * "." * string(run_id))

	# whether dry run
	if dry_run
		@goto dry_run_start
	end

	verbose_level = parse_verbose(verbose)

	if verbose_level == :all
		if p.info_before == "auto" || p.info_before == ""
			@info timestamp() * "Started: $(p.name)" run_id inputs outputs
		else
			@info timestamp() * p.info_before run_id inputs outputs
		end
	elseif verbose_level == :min
		if p.info_before == "auto" || p.info_before == ""
			@info timestamp() * "Started: $(p.name) [$(run_id)]"
		else
			@info timestamp() * p.info_before * " [$(run_id)]"
		end
	end

	# skip when done
	if skip_when_done && isfile(run_id_file) && isok(p.validate_outputs(outputs))
		if verbose_level == :all
			@warn timestamp() * "Skipped finished program: $(p.name)" run_id inputs outputs
		else
			@warn timestamp() * "Skipped finished program: $(p.name) [$(run_id)]"
		end
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
		@error timestamp() * "ProgramInputValidationError: $(p.name): fail to validate inputs (before running the main command)." validation_function=p.validate_inputs run_id inputs outputs
		rethrow(e)
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
		@error timestamp() * "ProgramPrerequisitesError: $(p.name): fail to run the prerequisites function (before running the main command)." prerequisites=p.prerequisites run_id inputs outputs
		rethrow(e)
		return false, outputs
	end

	# run the main command
	outputs_main = try
		p.main(inputs, outputs)
	catch e
		@error timestamp() * "ProgramRunningError: $(p.name): fail to run the main command." run_id inputs outputs
		rethrow(e)
		return false, outputs
	end

	# check type of outputs
	if !isa(outputs_main, Dict)
		@warn timestamp() * "ProgramMainReturnValue: $(p.name): the returned value of the main function is not a Dict, use the inferred output instead" run_id inputs outputs returned_outputs=outputs_main

		# check keys in outputs::Dict{String}
		check_outputs_keywords(p, outputs)
	else
		outputs = try
			check_outputs_keywords(p, outputs_main)
			outputs_main
		catch
			@warn timestamp() * "ProgramMainReturnValue: $(p.name): the returned Dict of the main function does not pass the keyword checks, use the inferred output instead" run_id inputs outputs returned_outputs=outputs_main
			outputs
		end
	end


	# confirmation: validate outputs
	try
		isok(p.validate_outputs(outputs)) || error("ProgramOutputValidationError: $(p.name): the validation function returns false.")
	catch e
		@error timestamp() * "ProgramOutputValidationError: $(p.name): fail to validate outputs (after running the main command)." validation_function=p.validate_outputs run_id inputs outputs
		rethrow(e)
		return false, outputs
	end

	try
		isok(p.wrap_up(inputs, outputs)) || error("ProgramWrapUpError: $(p.name): the wrap_up function returns false.")
	catch e
		@error timestamp() * "ProgramWrapUpError: $(p.name): fail to run the wrap_up function." wrap_up=p.wrap_up run_id inputs outputs
		rethrow(e)
		return false, outputs
	end

	# touch the run_id_file
	touch_run_id_file && touch(run_id_file)

	if verbose_level == :all
		if p.info_after == "auto" || p.info_after == ""
			@info timestamp() * "Finished: $(p.name)" run_id inputs outputs
		else
			@info timestamp() * p.info_after run_id inputs outputs
		end
	elseif verbose_level == :min
		if p.info_after == "auto" || p.info_after == ""
			@info timestamp() * "Finished: $(p.name) [$run_id]"
		else
			@info timestamp() * p.info_after * " [$run_id]"
		end
	end
	return true, outputs
end
