

## logging with dates
# JL Documentation page rendering may not compatible with LogginExtras, so use timestamp instead.
timestamp() = Dates.format(now(), "yyyy-mm-dd HH:MM:SS ")

## UUID generation
const UUID4 = uuid4(UUIDs.MersenneTwister(39871))

## simple functions
"""
    do_nothing() = nothing
    do_nothing(x) = nothing
    do_nothing(x, y) = nothing
"""
do_nothing() = nothing  # COV_EXCL_LINE
do_nothing(x) = nothing
do_nothing(x, y) = nothing

"""
    isok(x::Nothing) = true
    isok(x::Bool) = x
    isok(x::AbstractString) = true unless x is "" / n / no / null / f / false / 0
    isok(x::Any) = true  # default is true
"""
isok(x::Nothing) = true
isok(x::Bool) = x
isok(x::AbstractString) = !isempty(x) && !occursin(r"^n(o|ull)?$|^f(alse)?$|^0$"i, x)
isok(x) = true  # default is true


"""
    to_xxput_dict(p::Pair{String, V}) where V
    to_xxput_dict(p::Pair)
    to_xxput_dict(v::Vector{V}) where V <: Pair
    to_xxput_dict(d::Dict)

Convert inputs/outputs to Dict{String} in `run(p, inputs, outputs)`
"""
function to_xxput_dict(p::Pair{String, V}) where V
    Dict(p.first => p.second)
end
function to_xxput_dict(p::Pair)
    Dict(str(p.first) => p.second)
end
function to_xxput_dict(v::Vector{V}) where V <: Pair
    res = Dict{String,Any}()
    for p in v
        res[str(p.first)] = p.second
    end
    res
end
function to_xxput_dict(d::Dict)
    res = Dict{String,Any}()
    for p in d
        res[str(p.first)] = p.second
    end
    res
end
to_xxput_dict(d::Dict{String}) = d
to_xxput_dict(any) = throw(ErrorException("TypeError: cannot run Program: cannot convert to Dict{String}: inputs, outputs, or returned value of infer_outputs(inputs)"))


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
str(x::Regex) = x.pattern
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
to_cmd(any) = Cmd([str(any)])

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
