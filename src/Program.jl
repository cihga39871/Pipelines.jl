
"""
# Summary

	abstract type Program <: Any

# Subtypes

	CmdProgram
	JuliaProgram
"""
abstract type Program end

"""
	infer_outputs(p::Program, inputs::Dict{String})

Infer the default outputs from `p::Program` and `inputs::Dict{String}`.
"""
function infer_outputs(p::Program, inputs::Dict{String})
	p.infer_outputs(inputs)
end

function check_keywords(p::Program, inputs::Dict{String}, outputs::Dict{String})
	Set(p.inputs) == Set(keys(inputs)) || error("ProgramInputError: $(p.name): inputs provided and demanded not identical.")
	Set(p.outputs) == Set(keys(outputs)) || error("ProgramOutputError: $(p.name): outputs provided and demanded not identical.")
end
function check_outputs_keywords(p::Program, outputs::Dict{String})
	Set(p.outputs) == Set(keys(outputs)) || error("ProgramOutputError: $(p.name): outputs provided and demanded not identical.")
end

"""
	check_dependency(p::Program)

Check dependencies listed in `p.cmd_dependencies`.
"""
function check_dependency(p::Program)
	foreach(check_dependency, p.cmd_dependencies)
end

function generate_run_uuid(inputs::Dict{String}, outputs::Dict{String})
	all_cmd = ``
	foreach(x -> append!(all_cmd.exec, to_cmd(x).exec), values(inputs))
	foreach(x -> append!(all_cmd.exec, to_cmd(x).exec), values(outputs))
	uuid5(UUID4, string(all_cmd))
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
			check_data_type(inputs[keyword], p.input_types[i])
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
			check_data_type(outputs[keyword], p.output_types[i])
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
	xxputs_completion_and_check(p::Program, inputs::Dict{String}, outputs::Dict{String})

1. Check and complete `inputs` using types and values stored in `p`.

2. If `outputs` is empty, run `p.infer_outputs`.

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
		if isempty(outputs) && p.infer_outputs !== do_nothing
			outputs = p.infer_outputs(inputs)
		end
		outputs = outputs_completion(p::Program, outputs::Dict{String})
	end

	# check keyword consistency (keys in inputs and outputs compatible with the function)
	check_keywords(p, inputs, outputs)

	# parse <keyword> in String
	inputs, outputs = keyword_interpolation(inputs, outputs)
end

function Base.run(p::Program, inputs::Dict{String}, outputs::Dict{String}; kwarg...)
	run(p; inputs=inputs, outputs=outputs, kwarg...)
end

function Base.run(p::Program, inputs::Dict{String}; kwarg...)
	run(p; inputs=inputs, kwarg...)
	# if p.infer_outputs === do_nothing
	# 	if isempty(p.outputs)
	# 		run(p; inputs=inputs, kwarg...)
	# 	else
	# 		error("Cannot run Program '$(p.name)' without specifying outputs because it requires outputs but does not have a pre-defined `p.infer_outputs` function.")
	# 	end
	# else
	# 	outputs = infer_outputs(p, inputs)
	# 	run(p; inputs=inputs, outputs=outputs, kwarg...)
	# end
end
