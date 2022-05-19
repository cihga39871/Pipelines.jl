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
		name::String                            = "Command Program",
		id_file::String                         = "",
		info_before::String                     = "auto",
		info_after::String                      = "auto",
		cmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),
		inputs                                  = Vector{String}(),
		validate_inputs::Expr                   = do_nothing,  # vars of inputs
		infer_outputs::Expr                     = do_nothing,  # vars of inputs
		prerequisites::Expr                     = do_nothing,  # vars of inputs and outputs
		main::Expr                              = do_nothing,  # vars of inputs and outputs
		outputs                                 = Vector{String}(),
		validate_outputs::Expr                  = do_nothing,  # vars of outputs
		wrap_up::Expr                           = do_nothing   # vars of inputs and outputs
	) -> JuliaProgram

Julia program template. To run a `JuliaProgram`, use `run(::JuliaProgram; kwargs...).`

# Arguments

- `name::String`: Program name.

- `id_file::String`: The prefix of *run ID file*. To prevent from running the program with the same inputs and outputs twice, it will generate a unique *run ID file* after a successful run.

- `info_before::String`: Print it when the program is started.

- `info_after::String`: Print it when the program is finished.

- `cmd_dependencies::Vector{CmdDependency}`: Any command dependencies used in the program.

- `inputs` and `outputs`: Elements (or vectors containing elements) in the following format: (1) `keyword` (2) `keyword => data_type` (3) `keyword => default_value` (4) `keyword => default_value => data_type`.

  > `keyword` is an argument name, needs to be a `String`.
  >
  > `default_value` is optional. If set, users may not provide this argument when running. Elsewise, users have to provide it. Caution: `nothing` is preserved and means default value not set. If `String`, it can contain other keywords, but need to quote using '<>', such as `"<arg>.txt"`
  >
  > `data_type` is optional. If set, the value provided have to be this data type, or an error will throw.

  *HOW DOES THIS WORK?*

  > `JuliaProgram` stores a Julia function, with two arguments `inputs::Dict{String}, outputs::Dict{String}`. The keys of the two arguments should be set in `inputs::Vector{String}` and `outputs::Vector{String}`.
  > Caution: the returned value of `p.main` will be assigned to new `outputs` when the returned value is `Dict{String}` with proper keys. Otherwise, a warning will show and the returned value will be the original `outputs::Dict`.

- `validate_inputs::Expr`: A quoted code to validate inputs. Elements in `inputs` can be directly used as variables. If validation fail, throw error or return false. See details in [`quote_expr`](@ref)

- `infer_outputs::Expr`: A quoted code to infer outputs from inputs. Elements in `inputs` can be directly used as variables. Has to return a `Dict{String}("OUTPUT_VAR" => value)`. See details in [`quote_expr`](@ref)

- `prerequisites::Expr`: A quoted code to run just before the main command. It prepares necessary things, such as creating directories. Elements in `inputs` and `outputs` can be directly used as variables. See details in [`quote_expr`](@ref)

- `main::Expr`: The main julia code. Elements in `inputs` and `outputs` can be directly used as variables. See details in [`quote_expr`](@ref)

  > Caution: the returned value of `p.main` will be assigned to new `outputs`. Please ensure the variables in outputs are defined correctly, since it will return `outputs::Dict{String,Any}`.

- `validate_outputs::Expr`: A quoted code to validate outputs. Elements in `outputs` can be directly used as variables. If validation fail, throw error or return false. See details in [`quote_expr`](@ref)

- `wrap_up::Expr`: the last quoted code to run. Elements in `inputs` and `outputs` can be directly used as variables. See details in [`quote_expr`](@ref)

> **Compatibility of Pipelines < v0.8**:
>
> You can still pass `Function` to variables require `Expr`, but you **cannot** use the 'elements as variables' feature. The function should take `inputs::Dict{String}` and/or `outputs::Dict{String}` as variables, and you have to use traditional `inputs["VARNAME"]` in functions.
>
> From Pipelines v0.8, all `Expr` provided will be converted to `Function` automatically.

> **Debug: variable not found**
>
> Please refer to [`quote_expr`](@ref), section 'quote variables in other scopes.'

# Example

	p = JuliaProgram(
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

	# running the program: keyword arguments include keys of inputs and outputs
	success, new_out = run(p; a = `in1`, b = 2, c = "out", touch_run_id_file = false)

	@assert new_out != infer_outputs(p; a = `in1`, b = 2, c = "out")  # outputs will change to the returned value of main function, if the returned value is a Dict and pass `p.validate_outputs`

	# an old way to `run` program: need to create Dicts of inputs and outputs first.
	inputs = Dict("a" => `in1`, "b" => 2)
	outputs = "c" => "out"
	success, new_out = run(p, inputs, outputs; touch_run_id_file = false)
"""
function JuliaProgram(;
	name::String                            = "Julia Program",
	id_file::String                         = "",
	info_before::String                     = "auto",
	info_after::String                      = "auto",
	cmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),
	inputs                                  = Vector{String}(),
	validate_inputs::Union{Function,Expr}   = do_nothing,  # vars of inputs
	infer_outputs::Union{Function,Expr}     = do_nothing,  # vars of inputs
	prerequisites::Union{Function,Expr}     = do_nothing,  # vars of inputs and outputs
	main::Union{Function,Expr}              = do_nothing,  # vars of inputs and outputs
	outputs                                 = Vector{String}(),
	validate_outputs::Union{Function,Expr}  = do_nothing,  # vars of outputs
	wrap_up::Union{Function,Expr}           = do_nothing   # vars of inputs and outputs
)
	inputs, input_types, default_inputs = parse_default(inputs)
	outputs, output_types, default_outputs = parse_default(outputs)

	validate_inputs = quote_function(validate_inputs, inputs)
	infer_outputs = quote_function(infer_outputs, inputs)
	prerequisites = quote_function(prerequisites, inputs, outputs)
	validate_outputs = quote_function(validate_outputs, outputs)
	wrap_up = quote_function(wrap_up, inputs, outputs)

	main = quote_function(main, inputs, outputs; specific_return = :(outputs))

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
