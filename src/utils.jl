

## logging with dates
# JL Documentation page rendering may not compatible with LogginExtras, so use timestamp instead.
timestamp() = Dates.format(now(), "yyyy-mm-dd HH:MM:SS ")

## UUID generation
const UUID4 = uuid4(UUIDs.MersenneTwister(39871))

## simple functions
do_nothing() = nothing
do_nothing(x) = nothing
do_nothing(x, y) = nothing

isok(x::Nothing) = true
isok(x::Bool) = x
isok(x::AbstractString) = occursin(r"y(es)?|ok?|t(rue)?|^1$"i, x)
isok(x) = true  # default is true

## parse default inputs/outputs
"""
    parse_default(vec)

Parsing `inputs` and `outputs` when creating `Program` objects.

Return `xxputs::Vector{String}, xxput_types::Vector{DataType}, default_xxputs::Vector`.

## Valid `vec` element types

- `name::String`: no default value.

- `name::String => value`: set default value, except `value` is `nothing` or `DataType`.

- `name::String => value_type::DataType`: no default value, but value type.

- `name::String => value => value_type::DataType`: set default value and value type.

- `name::String => value_type::DataType => value`: set default value and value type.
"""
function parse_default(vec::Vector{String})
    xxput_types = DataType[Any for i = 1:length(vec)]
    default_xxputs = Nothing[nothing for i = 1:length(vec)]
    return vec, xxput_types, default_xxputs
end

function parse_default_element(ele::String)
    return ele, Any, nothing
end

function parse_default_element(ele::Pair{String,DataType})
    return ele.first, ele.second, nothing
end

function parse_default_element(ele::Pair{String,T}) where T
    return ele.first, Any, ele.second
end

function parse_default_element(ele::Pair{String,Pair{T,DataType}}) where T
    check_data_type(ele.second.first, ele.second.second)
    return ele.first, ele.second.second, ele.second.first
end

function parse_default_element(ele::Pair{String,Pair{DataType,T}}) where T
    check_data_type(ele.second.second, ele.second.first)
    return ele.first, ele.second.first, ele.second.second
end

function parse_default_element(ele::Pair{String,Pair{DataType,DataType}})
    if isa(ele.second.first, ele.second.second)
        return ele.first, ele.second.second, ele.second.first
    elseif isa(ele.second.second, ele.second.first)
        return ele.first, ele.second.first, ele.second.second
    else
        throw(ErrorException("DataTypeError: $(ele.second.first) and $(ele.second.second) are not inclusive."))
    end
end

function parse_default_element(ele::Any)
    throw(ErrorException("DataTypeError: $ele is not valid for a Program argument."))
end
function check_data_type(value, data_type::DataType)
    isa(value, data_type) || throw(ErrorException("DataTypeError: $value is not a $data_type type."))
end
## String/Cmd conversion
"""
    to_str(x) -> String
    str(x) -> String

Convert `x` to `String`.

- `x::Cmd`: remove backticks (return `string(x)[2:end-1]`).

- `x::Nothing`: return `""`.

- `x::Vector`: join elements with `"_"` as delim.

- `x::Any`: return `string(x)`.
"""
str(x::Cmd) = string(x)[2:end-1]
str(::Nothing) = ""
str(x::Vector) = join(str.(x), "_")
str(any) = string(any)
to_str = str

"""
    to_cmd(x) -> Cmd

Convert `x` to Cmd.

Exception: when `x::Nothing`, return `nothing::Nothing`.
"""
to_cmd(c::Cmd) = c
to_cmd(s::AbstractString) = Cmd([str(s)])
to_cmd(v::Vector{String}) = Cmd(v)
to_cmd(v::Vector{T}) where T = Cmd(str.(v))
to_cmd(::Nothing) = nothing
to_cmd(i::Number) = Cmd([str(i)])

ValidInputTypes = Union{Cmd, AbstractString, Vector, Number}

"""
    split(c::Cmd) = c.exec :: Vector{String}

Return splitted arguments of `Cmd`.
"""
function Base.split(c::Cmd)
    c.exec
end

## path functions
"""
    replaceext(path, replacement::AbstractString)

If the last component of a path contains a dot, leave everything before the dot as usual, and everything after the dot is replaced by `replacement`. Otherwise, `replacement` will be appended to `path`.

If `replacement` is empty, the last dot will be removed.
"""
function replaceext(path::String, replacement::AbstractString)
    a, b = splitext(path)
    if isempty(replacement)
        a
    elseif replacement[1] == '.'
        a * replacement
    else
        a * "." * replacement
    end
end
replaceext(path, replacement::AbstractString) = replaceext(to_str(path), replacement::String)

"""
    removeext(path)

If the last component of a path contains a dot, leave everything before the dot as usual, and everything including and after the dot is discarded.
"""
removeext(path::String) = splitext(path)[1]
removeext(path) = removeext(to_str(path))
