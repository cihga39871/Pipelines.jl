

function Base.display(p::CmdDependency)
	println("CmdDependency:")
	println("  exec             → ", p.exec)
	println("  test_args        → ", p.test_args)
	println("  validate_success → ", p.validate_success)
	println("  validate_stdout  → ", p.validate_stdout)
	println("  validate_stderr  → ", p.validate_stderr)
	println("  exit_when_fail   → ", p.exit_when_fail)
end

function Base.print(io::IO, p::CmdDependency)
	print(io, p.exec)
end

function Base.show(io::IO, p::CmdDependency)
	show(io, p.exec)
end

function Base.display(a::Arg)
	println("Arg:")
	println("  name        → ", a.name)
	println("  type        → ", a.type)
	println("  default     → ", a.default)
	println("  required    → ", a.required)
	println("  independent → ", a.independent)
end

function Base.display(args::Vector{Arg}; indent::Int = 2, summary::Bool = true)
	n_arg = length(args)
	if n_arg == 0
		if summary
			println(n_arg, "-element Vector{Arg}.")
		else
			println("<Arg empty>")
		end
		return
	end
	n_name = maximum([length(a.name) for a in args])
	n_type = maximum([length(string(a.type)) for a in args])
	summary && println(length(args), "-element Vector{Arg}:")
	indent_str = " " ^ indent
	for (i, a) in enumerate(args)
		space = (!summary && i == 1) ? "" : indent_str
		if a.required
			if a.independent
				s = "$space%-$(n_name)s :: %-$(n_type)s (required, independent)"
			else
				s = "$space%-$(n_name)s :: %-$(n_type)s (required)"
			end
			sout = @eval @sprintf($s, $(a.name), $(a.type))
		else
			if a.independent
				s = "$space%-$(n_name)s :: %-$(n_type)s (default = %s, independent)"
			else
				s = "$space%-$(n_name)s :: %-$(n_type)s (default = %s)"
			end
			sout = @eval @sprintf($s, $(a.name), $(a.type), $(string(a.default)))
		end
		println(sout)
	end
end

@eval function Base.display(p::CmdProgram)
    fs = $(fieldnames(CmdProgram))
    fs_string = $(map(string, fieldnames(CmdProgram)))
    max_byte = $(maximum(length, map(string, fieldnames(CmdProgram))))
    println("CmdProgram:")
    for (i,f) in enumerate(fs)
		print("  ", f, " " ^ (max_byte - length(fs_string[i])), " → ")
		if f === :arg_inputs || f === :arg_outputs
        	display(getfield(p, f); indent = max_byte + 5, summary = false)
		else
			println(getfield(p, f))
		end
    end
end

@eval function Base.display(p::JuliaProgram)
    fs = $(fieldnames(JuliaProgram))
    fs_string = $(map(string, fieldnames(JuliaProgram)))
    max_byte = $(maximum(length, map(string, fieldnames(JuliaProgram))))
    println("JuliaProgram:")
	for (i,f) in enumerate(fs)
		print("  ", f, " " ^ (max_byte - length(fs_string[i])), " → ")
		if f === :arg_inputs || f === :arg_outputs
        	display(getfield(p, f); indent = max_byte + 5, summary = false)
		else
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
