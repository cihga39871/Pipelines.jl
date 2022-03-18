

function Base.display(p::CmdDependency)
	print("CmdDependency\n  exec             →")
	display(p.exec)
	print("  test_args        →")
	display(p.test_args)
	println("  validate_success → $(p.validate_success)\n  validate_stdout  → $(p.validate_stdout)\n  validate_stderr  → $(p.validate_stderr)\n  exit_when_fail   → $(p.exit_when_fail)")
end

function Base.print(io::IO, p::CmdDependency)
	print(io, p.exec)
end

function Base.show(io::IO, p::CmdDependency)
	show(io, p.exec)
end

@eval function Base.display(p::CmdProgram)
    fs = $(fieldnames(CmdProgram))
    fs_string = $(map(string, fieldnames(CmdProgram)))
    max_byte = $(maximum(length, map(string, fieldnames(CmdProgram))))
    println("$CmdProgram:")
    for (i,f) in enumerate(fs)
		if f == :inputs
			print("  ", f, " " ^ (max_byte - length(fs_string[i])), " → ")
			display_xxputs(max_byte, p.inputs, p.input_types, p.default_inputs)
		elseif f == :outputs
			print("  ", f, " " ^ (max_byte - length(fs_string[i])), " → ")
			display_xxputs(max_byte, p.outputs, p.output_types, p.default_outputs)
		elseif f in Symbol[:input_types, :default_inputs, :output_types, :default_outputs]
		else
			print("  ", f, " " ^ (max_byte - length(fs_string[i])), " → ")
	        println(getfield(p, f))
		end
    end
end

@eval function Base.display(p::JuliaProgram)
    fs = $(fieldnames(JuliaProgram))
    fs_string = $(map(string, fieldnames(JuliaProgram)))
    max_byte = $(maximum(length, map(string, fieldnames(JuliaProgram))))
    println("$JuliaProgram:")
    for (i,f) in enumerate(fs)
		if f == :inputs
			print("  ", f, " " ^ (max_byte - length(fs_string[i])), " → ")
			display_xxputs(max_byte, p.inputs, p.input_types, p.default_inputs)
		elseif f == :outputs
			print("  ", f, " " ^ (max_byte - length(fs_string[i])), " → ")
			display_xxputs(max_byte, p.outputs, p.output_types, p.default_outputs)
		elseif f in Symbol[:input_types, :default_inputs, :output_types, :default_outputs]
		else
			print("  ", f, " " ^ (max_byte - length(fs_string[i])), " → ")
	        println(getfield(p, f))
		end
    end
end

function display_xxputs(max_bype::Int, xxputs::Vector{String}, xxput_types::Vector{Type}, default_xxputs::Vector)
	n = length(xxputs)
	if n == 0
		println("<empty>")
		return
	end
	max_bype_xxputs = maximum(length, xxputs)
	max_bype_xxput_types = maximum(length, map(string, xxput_types))
	for i in 1:n
		name = xxputs[i]
		type = string(xxput_types[i])
		default = default_xxputs[i]
		println(i == 1 ? "" : " " ^ (max_bype + 5),
			"\"$(name)\" ", " " ^ (max_bype_xxputs - length(name)),
			":: $(type) ", " " ^ (max_bype_xxput_types - length(type)),
			isnothing(default) ? "(required)" : "(default: $default)"
		)
	end
	return
end
