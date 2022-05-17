
"""
# Summary

	abstract type Program <: Any

# Subtypes

	CmdProgram
	JuliaProgram
"""
abstract type Program end

"""
	infer_outputs(p::Program, inputs)
	infer_outputs(p::Program, inputs, outputs)

Infer the default outputs from `p::Program` and `inputs::Dict{String}`.
"""
function infer_outputs(p::Program, inputs, outputs = Dict())
	i, o = xxputs_completion_and_check(p, inputs, outputs)
	o
end

function check_keywords(p::Program, inputs::Dict{String}, outputs::Dict{String})
	Set(p.inputs) == Set(keys(inputs)) || error("ProgramInputError: $(p.name): inputs provided and demanded not identical.")
	Set(p.outputs) == Set(keys(outputs)) || error("ProgramOutputError: $(p.name): outputs provided and demanded not identical.")
end
function check_outputs_keywords(p::Program, outputs::Dict{String})
	Set(p.outputs) == Set(keys(outputs)) || error("ProgramOutputError: $(p.name): outputs provided and demanded not identical.")
end

"""
	check_dependency(p::Program; exit_when_fail::Bool=true)

Check dependencies listed in `p.cmd_dependencies`.
"""
function check_dependency(p::Program; exit_when_fail::Bool=true)
	length(p.cmd_dependencies) == 0 && (return true)
	final_res = true
	for dep in p.cmd_dependencies
		res = check_dependency(dep; exit_when_fail = exit_when_fail)
		if !res
			final_res = false
		end
	end
	final_res
end

"""
    check_dependency(m::Module = @__MODULE__; exit_when_fail = true, verbose = true)

Check all `CmdDependency` and `Program` under `m::Module`.
"""
function check_dependency(m::Module = @__MODULE__; exit_when_fail = true, verbose = true)
	verbose && (@info "Dependencies Status:")

    deps = filter!(x -> isdefined(m, x) && (getfield(m, x) isa CmdDependency || getfield(m, x) isa Program), names(m))
	length(deps) == 0 && (return true)

	final_res = true
	max_len = maximum([length(String(d)) for d in deps])
	for dep_name in deps
		nspace = max_len - length(String(dep_name))
		space = " " ^ nspace
		dep = getfield(m, dep_name)
		res = check_dependency(dep; exit_when_fail = exit_when_fail)
		status = res ? "  OK" : "FAIL"
		if verbose
			if dep isa CmdDependency
				@info "  $status  $m.$dep_name  $space$(dep.exec)"
			else  # Program
				@info "  $status  $m.$dep_name"
			end
		end

		if !res
			final_res = false
		end
	end
	return final_res
end

"""
    status_dependency(m::Module = @__MODULE__; exit_when_fail = false, verbose = true)

Check all `CmdDependency` and `Program` under `m::Module`. Similar to `check_dependency`, but do not `exit_when_fail` by default.
"""
function status_dependency(m::Module = @__MODULE__; exit_when_fail = false, verbose = true)
	check_dependency(m; exit_when_fail = exit_when_fail, verbose = verbose)
end

function generate_run_uuid(inputs::Dict{String}, outputs::Dict{String})
	out_uuid = UUID4
	for d in (inputs, outputs)
		sd = sort(OrderedDict(d))
		for (k,v) in sd
			out_uuid = uuid5(out_uuid, string(k, ":", v))
		end
	end
	return out_uuid
end

"""
	inputs_completion(p::Program, inputs::Dict{String})

complete missing keywords in `inputs` and `outputs`.
"""
function inputs_completion(p::Program, inputs::Dict{String})
	value_type = fieldtypes(eltype(inputs))[2]  # value type of inputs
	for (i, keyword) in enumerate(p.inputs)
		if haskey(inputs, keyword)
			# check types
			value = convert_data_type(inputs[keyword], p.input_types[i])
			if !(value isa value_type)
				inputs = convert(Dict{String,Any}, inputs)
			end  # it is ok to replace in for-loop
			inputs[keyword] = value
		else
			# not provided, check default
			default = p.default_inputs[i]
			if isnothing(default)
				# no default
				throw(ErrorException("ArgumentError: Program '$(p.name)' requires '$keyword' in inputs, but it is not provided."))
			else
				if !(default isa value_type)
					inputs = convert(Dict{String,Any}, inputs)
				end  # it is ok to replace in for-loop
				inputs[keyword] = default
			end
		end
	end
	inputs
end
function outputs_completion(p::Program, outputs::Dict{String})
	value_type = fieldtypes(eltype(outputs))[2]  # value type of outputs
	for (i, keyword) in enumerate(p.outputs)
		if haskey(outputs, keyword)
			# check types
			value = convert_data_type(outputs[keyword], p.output_types[i])
			if !(value isa value_type)
				outputs = convert(Dict{String,Any}, outputs)
			end  # it is ok to replace in for-loop
			outputs[keyword] = value
		else
			# not provided, check default
			default = p.default_outputs[i]
			if isnothing(default)
				# no default
				throw(ErrorException("ArgumentError: Program '$(p.name)' requires $keyword in outputs, but it is not provided."))
			else
				if !(default isa value_type)
					outputs = convert(Dict{String,Any}, outputs)
				end  # it is ok to replace in for-loop
				outputs[keyword] = default
			end
		end
	end
	outputs
end

"""
	keyword_interpolation(inputs::Dict{String}, outputs::Dict{String})

Interpolate <keyword> in `String`.
"""
function keyword_interpolation(inputs::Dict{String}, outputs::Dict{String})
	allputs = Dict(union(inputs, outputs))
	keyword_interpolation(allputs::Dict{String})
	for key in keys(inputs)
		inputs[key] = allputs[key]
	end
	for key in keys(outputs)
		outputs[key] = allputs[key]
	end
	inputs, outputs
end
function keyword_interpolation(allputs::Dict{String})
	for key in keys(allputs)
		keyword_interpolation(allputs, key, 1)
	end
	allputs
end
function keyword_interpolation(allputs::Dict{String}, key::String, n_recursion::Int)
	value = allputs[key]
	keywords = find_keywords(value)
	isempty(keywords) && return allputs # no need to interpolate
	for keyword in keywords
		keyword_value = get(allputs, keyword, nothing)
		isnothing(keyword_value) && continue # keyword not match, ignore
		if isempty(find_keywords(keyword_value))
			allputs[key] = replace(allputs[key], "<$keyword>" => str(keyword_value))
		else
			n_recursion > 20 && throw(ErrorException("ProgramArgumentError: too many recursion occurs in keyword interpolation."))
			keyword_interpolation(allputs::Dict{String}, keyword::String, n_recursion + 1)
		end
	end
	allputs
end

function find_keywords(value::T) where T <: AbstractString
	m = eachmatch(r"<([^<>]+)>", value)
	String[i.captures[1] for i in m]
end
find_keywords(not_string) = []

"""
	xxputs_completion_and_check(p::Program, inputs, outputs)

1. Check and complete `inputs` using types and values stored in `p`.

2. Run `p.infer_outputs` if defined, and then merge it and outputs (user-input keys are kept).

3. Check and complete `outputs` using types and values stored in `p`.

4. Check keyword consistency using `p`.

5. Interpolate <keyword> in String in completed `inputs` and `outputs`.

6. Return inputs and outputs.
"""
function xxputs_completion_and_check(p::Program, inputs::Dict{String}, outputs::Dict{String})
	inputs = inputs_completion(p::Program, inputs::Dict{String})

	if isempty(p.outputs)
		# do nothing to outputs
	else
		if p.infer_outputs !== do_nothing
			outputs_from_infer = to_xxput_dict(p.infer_outputs(inputs))
			outputs = merge(outputs_from_infer, outputs)
		end
		outputs = outputs_completion(p::Program, outputs::Dict{String})
	end

	# check keyword consistency (keys in inputs and outputs compatible with the function)
	check_keywords(p, inputs, outputs)

	# parse <keyword> in String
	inputs, outputs = keyword_interpolation(inputs, outputs)
end

function xxputs_completion_and_check(p::Program, inputs, outputs)
	xxputs_completion_and_check(p, to_xxput_dict(inputs), to_xxput_dict(outputs))
end

function Base.run(p::Program, inputs, outputs; kwarg...)
	run(p; inputs=inputs, outputs=outputs, kwarg...)
end

function Base.run(p::Program, inputs; kwarg...)
	run(p; inputs=inputs, kwarg...)
end

function Base.run(p::Program;
	dir::AbstractString = "", retry::Int = 0,
	stdout = nothing, stderr = nothing, stdlog = stderr, append::Bool = false,
	kwarg...
)
	if dir != ""
		dir_backup = pwd()
		dir = abspath(dir)
		cd(dir) # go to working directory
	end

	n_try = 0
	local res
	while n_try <= retry
		res = redirect_to_files(stdout, stderr, stdlog; mode = append ? "a+" : "w+") do
			if n_try > 0
				@warn "Retry $(p.name) ($n_try/$retry)"
			end
			_run(p; dir = dir, kwarg...)
		end
		res isa StackTraceVector || break  # res isa StackTraceVector means failed, need retry.
		n_try += 1
	end

	if dir != ""
		cd(dir_backup)
	end
	res
end

"""
	prog_run(p::Program; kwargs...)

Run Program (CmdProgram or JuliaProgram).

Return `(success::Bool, outputs::Dict{String})`

  > If `p isa JuliaProgram`, `outputs` **will be overwritten by the returned value of** `p.main` ***only*** when the returned value is a `Dict{String}` and passes `p.validate_outputs`. See more at [`JuliaProgram`](@ref).

# Keyword Arguments:

- elements in `p.inputs` and `p.outputs`.

- `dir::AbstractString = ""`: working directory to run the program and store `run_id_file`.

- `check_dependencies::Bool = true`: check dependencies for `p` (`p.cmd_dependencies`).

- `skip_when_done::Bool = true`: Skip running the program and return `true` if it has been done before (the `run_id_file` exists and `p.validate_outputs(outputs)` passes.)

- `touch_run_id_file::Bool = true`: If `true`, touch a unique run ID file, which indicate the program is successfully run with given inputs and outputs. If `false`, the next time running the program, `skip_when_done=true` will not take effect.

- `verbose = true`: If `true` or `:all`, print all info and error messages. If `:min`, print minimum info and error messages. If `false` or `:none`, print error messages only.

- `retry::Int = 0`: If failed, retry for INT times.

- `dry_run::Bool = false`: do not run the program, return `(command::AbstractCmd, run_id_file::String)` for CmdProgram, or `(inferred_outputs::Dict{String}, run_id_file::String)` for JuliaProgram.

- `stdout`, `stderr`, `stdlog` and `append`: Redirect the program outputs to files. `stdlog` is the Julia logging of `@info`, `@warn`, `@error`, etc. Caution: If `p isa CmdProgram` and the original command (`p.cmd`) has redirection, arguments defined here might not be effective for the command.

!!! note
	Redirecting in Julia are not thread safe, so unexpected redirection might be happen if you are running programs in different `Tasks` or multi-thread mode.

### Workflow

1. Go to the working directory. Establish redirection. (`dir`, `stdout`, `stderr`, `stdlog`, `append`).

2. Validate compatibility between `p` and `inputs/outputs`.

3. Check whether the program has run before. (`skip_when_done`, `p.validate_outputs(outputs)`)

4. Check command dependencies. (`check_dependencies`, `p.cmd_dependencies`)

5. Validate `inputs`. (`p.validate_inputs(inputs)`)

6. [CmdProgram only] Generate runnable command from `p` and `inputs/outputs`. (`stdout`, `stderr`, `append`)

7. Preparing before running main command. (`p.prerequisites(inputs, outputs)`)

8. Run command [CmdProgram] or the main function [JuliaProgram].

9. If `p isa CmdProgram`, validate `outputs` only. If `p isa JuliaProgram`, validate the returned value of the main function. If pass, `outputs` will ***overwritten by the returned value***. Otherwise, the original `outputs` is ***kept***. (`p.validate_outputs(outputs)`)

10. Wrap up. (`p.wrap_up(inputs, outputs)`)

11. Success, touch run id file, and return `(success::Bool, outputs::Dict{String})`. (`touch_run_id_file::Bool`)


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

	# for CmdProgram, outputs are inferred before running the main command, however,
	# for JuliaProgram, outputs will change to the returned value of main function, if the returned value is a Dict and pass `p.validate_outputs`
	@assert new_out != outputs

	# an old way to `run` program: need to create Dicts of inputs and outputs first.
	inputs = Dict("a" => `in1`, "b" => 2)
	outputs = "c" => "out"
	success, new_out = run(p, inputs, outputs; touch_run_id_file = false)
"""
function prog_run(p::Program; args...)
    inputs, outputs, kws = parse_program_args(p::Program; args...)
	if isempty(inputs)

		if isempty(outputs)
			if length(kws) == 0
		        run(p)
		    else
		        run(p; kws...)
		    end
		else
			if length(kws) == 0
				run(p, inputs, outputs)
			else
				run(p, inputs, outputs; kws...)
			end
		end
	else
		if isempty(outputs)
			if length(kws) == 0
		        run(p, inputs)
		    else
		        run(p, inputs; kws...)
		    end
		else
			if length(kws) == 0
				run(p, inputs, outputs)
			else
				run(p, inputs, outputs; kws...)
			end
		end
	end
end

"""
	parse_program_args(p::Program; kwargs...)

Classify `args...` to inputs and outputs of `p` or other keyword arguments.

Return (inputs::Dict{String}, outputs::Dict{String}, other_kwargs::Tuple)
"""
function parse_program_args(p::Program; args...)
    inputs = Dict{String,Any}()
    outputs = Dict{String,Any}()
    kw_indices = Vector{Int}()
    for (i, arg) in enumerate(args)
        k, v = arg
        k_str = string(k)
        if k_str in p.inputs
            inputs[k_str] = v
        elseif k_str in p.outputs
            outputs[k_str] = v
        else
            push!(kw_indices, i)
        end
    end
    kws = collect(args)[kw_indices]
    return inputs, outputs, kws
end

function parse_verbose(verbose::Symbol)
	if verbose in (:all, :min, :none)
		return verbose
	elseif verbose == :minimum
		return :min
	elseif verbose in (:max, :maximum, :full, :yes, :true)
		return :all
	elseif verbose in (:no, :nothing, :null, :false)
		return :none
	else
		@error "Cannot determine verbose level ($verbose). Set to :all. Accepted verbose options are true, false, :all, :min, and :none."
		return :all
	end
end
parse_verbose(verbose::Bool) = verbose ? :all : :none
parse_verbose(::Nothing) = :none
parse_verbose(verbose::AbstractString) = parse_verbose(Symbol(verbose))
