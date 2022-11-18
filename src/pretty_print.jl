

function Base.show(io::IO, ::MIME"text/plain", p::CmdDependency)
    println(io, "CmdDependency:")
    println(io, "  exec             → ", p.exec)
    println(io, "  test_args        → ", p.test_args)
    println(io, "  validate_success → ", p.validate_success)
    println(io, "  validate_stdout  → ", p.validate_stdout)
    println(io, "  validate_stderr  → ", p.validate_stderr)
    println(io, "  exit_when_fail   → ", p.exit_when_fail)
end
function Base.show(io::IO, p::CmdDependency)
    show(io, p.exec)
end
function Base.print(io::IO, p::CmdDependency)
    print(io, p.exec)
end

function Base.show(io::IO, ::MIME"text/plain", a::Arg)
    println(io, "Arg:")
    println(io, "  name        → ", a.name)
    println(io, "  type        → ", a.type)
    println(io, "  default     → ", a.default)
    println(io, "  required    → ", a.required)
    println(io, "  independent → ", a.independent)
end
function Base.show(io::IO, a::Arg)
    type_str = a.type == Any ? "" : "::$(a.type)"
    if a.required
        if a.independent
            s = "$(a.name)$(type_str) (required, independent)"
        else
            s = "$(a.name)$(type_str) (required)"
        end
    else
        if a.independent
            s = "$(a.name)$(type_str) (independent, default = $(a.default))"
        else
            s = "$(a.name)$(type_str) (default = $(a.default))"
        end
    end
    print(io, s)
end

function Base.show(io::IO, ::MIME"text/plain", args::Vector{T}; indent::Int = 2, summary::Bool = true) where T <: Arg
    n_arg = length(args)
    if n_arg == 0
        if summary
            println(io, n_arg, "-element Vector{Arg}.")
        else
            println(io, "<Arg empty>")
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
                s = "$space%-$(n_name)s :: %-$(n_type)s (independent, default = %s)"
            else
                s = "$space%-$(n_name)s :: %-$(n_type)s (default = %s)"
            end
            sout = @eval @sprintf($s, $(a.name), $(a.type), $(string(a.default)))
        end
        println(io, sout)
    end
end
function show_vector(io::IO, ::MIME"text/plain", args::Vector{T}; indent::Int = 2, summary::Bool = true) where T
    n_arg = length(args)
    if n_arg == 0
        if summary
            println(io, n_arg, "-element Vector{$T}.")
        else
            println(io, "<empty>")
        end
        return
    end
    summary && println(length(args), "-element Vector{$T}:")
    indent_str = " " ^ indent
    for (i, a) in enumerate(args)
        space = (!summary && i == 1) ? "" : indent_str
        println(io, space, a)
    end
end

@eval function Base.show(io::IO, mime::MIME"text/plain", p::CmdProgram)
    fs = $(fieldnames(CmdProgram))
    fs_string = $(map(string, fieldnames(CmdProgram)))
    max_byte = $(maximum(length, map(string, fieldnames(CmdProgram))))
    println(io, "CmdProgram:")
    for (i,f) in enumerate(fs)
        print(io, "  ", f, " " ^ (max_byte - length(fs_string[i])), " → ")
        if f === :arg_inputs || f === :arg_outputs
            show(io, mime, getfield(p, f); indent = max_byte + 5, summary = false)
        elseif f in (:cmd_dependencies, :arg_forward)
            show_vector(io, mime, getfield(p, f); indent = max_byte + 5, summary = false)
        else
            println(io, getfield(p, f))
        end
    end
end
@eval function Base.show(io::IO, mime::MIME"text/plain", p::JuliaProgram)
    fs = $(fieldnames(JuliaProgram))
    fs_string = $(map(string, fieldnames(JuliaProgram)))
    max_byte = $(maximum(length, map(string, fieldnames(JuliaProgram))))
    println(io, "JuliaProgram:")
    for (i,f) in enumerate(fs)
        print(io, "  ", f, " " ^ (max_byte - length(fs_string[i])), " → ")
        if f === :arg_inputs || f === :arg_outputs
            show(io, mime, getfield(p, f); indent = max_byte + 5, summary = false)
        elseif f in (:cmd_dependencies, :arg_forward)
            show_vector(io, mime, getfield(p, f); indent = max_byte + 5, summary = false)
        else
            println(io, getfield(p, f))
        end
    end
end

function Base.show(io::IO, p::CmdProgram)
    print(io, "CmdProgram($(p.name), $(p.cmd))")
end
function Base.show(io::IO, p::JuliaProgram)
    print(io, "JuliaProgram($(p.name), inputs = $(p.inputs), outputs = $(p.outputs))")
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
