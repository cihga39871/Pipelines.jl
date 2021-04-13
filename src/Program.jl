
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


function Base.run(p::Program, inputs::Dict{String}, outputs::Dict{String}; kwarg...)
	run(p; inputs=inputs, outputs=outputs, kwarg...)
end

function Base.run(p::Program, inputs::Dict{String}; kwarg...)
	if p.infer_outputs === do_nothing
		if isempty(p.outputs)
			run(p; inputs=inputs, kwarg...)
		else
			error("Cannot run Program '$(p.name)' without specifying outputs because it requires outputs but does not have a pre-defined `p.infer_outputs` function.")
		end
	else
		outputs = infer_outputs(p, inputs)
		run(p; inputs=inputs, outputs=outputs, kwarg...)
	end
end
