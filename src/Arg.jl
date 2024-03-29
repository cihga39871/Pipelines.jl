
### Reserved keys
"""
    const RESERVED_KEY_SET = Set(["name", "user", "ncpu", "mem",
        "schedule_time", "wall_time", "priority", "dependency",
        "stdout", "stderr", "stdlog", "append", "dir", "inputs", "outputs",
        "check_dependencies", "skip_when_done", "touch_run_id_file",
        "verbose", "retry", "dry_run"])

Reserved keys that cannot be used in inputs and outputs.
"""
const RESERVED_KEY_SET = Set(["name", "user", "ncpu", "mem", "schedule_time", "wall_time", "priority", "dependency", "stdout", "stderr", "stdlog", "append", "dir", "inputs", "outputs", "check_dependencies", "skip_when_done", "touch_run_id_file", "verbose", "retry", "dry_run"])


### Arg

"""
```julia

Arg(name)
Arg(name => default)
Arg(name => type::Type)
Arg(name => default => type::Type)
Arg(name => type::Type => default)

Arg(name::Union{String,Symbol}, type::Type = Any, default = nothing;
    required::Bool = isnothing(default),
    independent::Bool = name isa Symbol
)

struct Arg{type,DefaultType}
    name::String
    type::Type
    default::DefaultType
    required::Bool
    independent::Bool
end
```

`Arg` stores the settings of inputs and outputs in `Program`.

- `name`: name of Arg.

- `type`: allowed type.

- `default`: default value.

- `required = isnothing(default)`: if true, the `Arg` has to be provided by users.

- `independent = isa(name, Symbol)`: if true, the argument does not change the results of a Program, such as "nthreads", "memory". Independent args have no effect on run id.


## Valid `pair` types

- `name`: no default value.

- `name => default`: set default value, except `value` is `nothing` (default value not set).

- `name => type`: no default value, but value type is restricted.

- `name => default => type`: set default value and value type.

- `name => type => default`: set default value and value type.


!!! tip "An edge situation"
    To create an argument with a default value of `nothing`, you cannot use `=>`. Instead, this works:
    ```julia
    p = JuliaProgram(
        inputs = [
            Arg("ARG_NAME", nothing; required = false),
            "OTHER_ARG" => String
        ]
    )
    ```
"""
struct Arg{AllowedType,DefaultType}
    name::String
    type::Type
    default::DefaultType
    required::Bool
    independent::Bool
    function Arg(name::String, type::Type, default::DefaultType, required::Bool, independent::Bool) where DefaultType
        if name in RESERVED_KEY_SET
            error("Program: cannot use the reserved name ($name) in inputs or outputs. Please use another name. All reserved names includes in the $RESERVED_KEY_SET")
        end
        default = convert_data_type(default, type)
        # if !(default isa type || default isa Nothing)
        #     throw(TypeError(:Arg, "checking default value and allowed type", type, default))
        # end
        new{type,DefaultType}(name, type, default, required, independent)
    end
end

function Arg(name::String, type::Type = Any, default = nothing; required::Bool = isnothing(default), independent::Bool = false)
    Arg(name, type, default, required, independent)
end
function Arg(name::AbstractString, type::Type = Any, default = nothing; required::Bool = isnothing(default), independent::Bool = false)
    Arg(String(name), type, default, required, independent)
end
function Arg(name::Symbol, type::Type = Any, default = nothing; required::Bool = isnothing(default), independent::Bool = true)
    Arg(String(name), type, default, required, independent)
end

function Arg(name::String, default; required::Bool = isnothing(default), independent::Bool = false)
    Arg(name, Any, default, required, independent)
end
function Arg(name::AbstractString, default; required::Bool = isnothing(default), independent::Bool = false)
    Arg(String(name), Any, default, required, independent)
end
function Arg(name::Symbol, default; required::Bool = isnothing(default), independent::Bool = true)
    Arg(String(name), Any, default, required, independent)
end


## search
function Base.in(name::AbstractString, args::Vector{Arg})
    for a in args
        if name == a.name
            return true
        end
    end
    return false
end

## parse default inputs/outputs
function throw_invalid_xxputs(any)
    throw(ErrorException("TypeError: Invalid value: $any. Elements of inputs and outputs can only be one of `name::Union{String,Symbol}`, `name => default_value`, `name => data_type`, `name => default_value => data_type`, `name => data_type => default_value`."))
end

"""
    parse_arg(v)

Parsing `inputs` and `outputs` when creating `Program` objects.

Return `Vector{Arg}`.

## Valid `v` element types

- `name`: no default value.

- `name => value`: set default value, except `value` is `nothing` (default value not set).

- `name => value_type::Type`: no default value, but value type.

- `name => value => value_type::Type`: set default value and value type.

- `name => value_type::Type => value`: set default value and value type.
"""
parse_arg(v::Vector{String}) = Arg[Arg(ele) for ele in v]
parse_arg(v::Vector) = Arg[Arg(ele) for ele in v]

parse_arg(s::String) = Arg[Arg(s)]
parse_arg(s::AbstractString) = Arg[Arg(s)]
parse_arg(p::Pair) = Arg[Arg(p)]
parse_arg(d::Dict) = Arg[Arg(p) for p in d]
parse_arg(s::Symbol) = Arg[Arg(s)]
parse_arg(a::Arg) = Arg[a]
parse_arg(any) = throw_invalid_xxputs(any)

function Arg(ele::Pair)
    name = ele.first
    type, value = parse_arg_info(ele.second)
    Arg(name, type, value)
end
Arg(a::Arg) = a
Arg(any) = throw_invalid_xxputs(any)


function parse_arg_info(data_type::Type)
    data_type, nothing
end
function parse_arg_info(p::Pair)
    parse_arg_info_pair(p.first, p.second)
end
function parse_arg_info(value)
    Any, value
end
function parse_arg_info_pair(a::Type, b::Type)
    if a isa b
        b, a
    elseif b isa a
        a, b
    else
        throw(ErrorException("TypeError: $(a) and $(b) are not inclusive."))
    end
end
parse_arg_info_pair(a, b::Type) = b, convert_data_type(a, b)
parse_arg_info_pair(a::Type, b) = a, convert_data_type(b, a)

convert_data_type(value::Nothing, data_type::Type) = nothing

function convert_data_type(value, data_type::Type)
    if isa(value, data_type)
        return value
    else
        try
            convert(data_type, value)
        catch
            throw(ErrorException("TypeError: cannot convert $value to $data_type."))
        end
    end
end

######## arg forward

"""
    FORWARD_KEY_SET = Set([:name, :user, :ncpu, :mem])

`Program` objects has a field `arg_forward`. It can forward args from inputs and outputs to `JobSchedulers.Job()`, only supporting keyword arguments in `Pipelines.FORWARD_KEY_SET`. `arg_forward` accepts elements (or vectors containing elements) in the following format: `\"arg_of_inputs_or_outputs\" => :key_in_FORWARD_KEY_SET`."
"""
const FORWARD_KEY_SET = Set([:name, :user, :ncpu, :mem])

function throw_invalid_type_of_arg_forward(any)
    throw(ErrorException("TypeError: Invalid value: $any. `arg_forward` can forward args from inputs and outputs to `JobSchedulers.Job()`, only supporting keyword arguments in `Pipelines.FORWARD_KEY_SET`: `$FORWARD_KEY_SET`. `arg_forward` accepts elements (or vectors containing elements) in the following format: `\"arg_of_inputs_or_outputs\" => :key_in_FORWARD_KEY_SET`."))
end

parse_arg_forward(v::Vector{Pair{String,Symbol}}) = v
parse_arg_forward(v::Vector) = Pair{String,Symbol}[parse_arg_forward_element(p) for p in v]
parse_arg_forward(p::Pair) = [parse_arg_forward_element(p)]
parse_arg_forward(d::Dict) = parse_arg_forward(collect(d))
parse_arg_forward(any) = throw_invalid_type_of_arg_forward(any)

parse_arg_forward_element(p::Pair{String,Symbol}) = p
parse_arg_forward_element(p::Pair) = String(p.first) => Symbol(p.second)
parse_arg_forward_element(any) = throw_invalid_type_of_arg_forward(any)

function namein(x::String, arg_inputs::Vector{Arg}, arg_outputs::Vector{Arg})
    for arg in arg_inputs
        if arg.name == x
            return true
        end
    end
    for arg in arg_outputs
        if arg.name == x
            return true
        end
    end
    false
end

function check_arg_forward(arg_forward::Vector{Pair{String,Symbol}},arg_inputs::Vector{Arg}, arg_outputs::Vector{Arg})
    isempty(arg_forward) && (return true)
    for p in arg_forward
        if p.second in FORWARD_KEY_SET
            if !namein(p.first, arg_inputs, arg_outputs)
                throw(ArgumentError("invalid arg_forward element \"$(p.first)\" => :$(p.second): \"$(p.first)\" is not found in inputs or outputs."))
            end
        else
            throw(ArgumentError("invalid arg_forward element \"$(p.first)\" => :$(p.second): :$(p.second) cannot be forwarded to other functions. Only names in $FORWARD_KEY_SET can be forwarded."))
        end
    end
    true
end