

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
function throw_invalid_xxputs(any)
	throw(ErrorException("DataTypeError: Elements of inputs and outputs can only be one of `name::Union{String,Symbol}`, `name => default_value`, `name => data_type`, `name => default_value => data_type`, `name => data_type => default_value`. Invalid value: $any"))
end

"""
    parse_default(v)

Parsing `inputs` and `outputs` when creating `Program` objects.

Return `xxputs::Vector{String}, xxput_types::Vector{DataType}, default_xxputs::Vector`.

## Valid `v` element types

- `name::String`: no default value.

- `name::String => value`: set default value, except `value` is `nothing` (default value not set).

- `name::String => value_type::DataType`: no default value, but value type.

- `name::String => value => value_type::DataType`: set default value and value type.

- `name::String => value_type::DataType => value`: set default value and value type.
"""
function parse_default(v::Vector{String})
    xxput_types = DataType[Any for i = 1:length(v)]
    default_xxputs = Nothing[nothing for i = 1:length(v)]
    return v, xxput_types, default_xxputs
end

function parse_default(v::Vector)
    n = length(v)
    xxputs = Vector{String}(undef, n)
    xxput_types = Vector{DataType}(undef, n)
    default_xxputs = Vector{Any}(undef, n)
    for (i, ele) in enumerate(v)
        name, type, default = parse_default_element(ele)
        xxputs[i] = name
        xxput_types[i] = type
        default_xxputs[i] = default
    end
    return xxputs, xxput_types, default_xxputs
end

parse_default(s::String) = parse_default([s])
parse_default(s::AbstractString) = parse_default([str(s)])
parse_default(p::Pair) = parse_default([p])
parse_default(d::Dict) = parse_default([p for p in d])
parse_default(s::Symbol) = parse_default([str(s)])
parse_default(any) = throw_invalid_xxputs(any)

function parse_default_element(ele::String)
    return ele, Any, nothing
end
function parse_default_element(ele::AbstractString)
    return str(ele), Any, nothing
end
function parse_default_element(ele::Symbol)
    return str(ele), Any, nothing
end
function parse_default_element(ele::Pair)
	if ele.first isa AbstractString || ele.first isa Symbol
	    name = str(ele.first)
		data_type, value = parse_arg_info(ele.second)
		name, data_type, value
	else
		throw_invalid_xxputs(ele)
	end
end
parse_default_element(any) = throw_invalid_xxputs(any)

function parse_arg_info(data_type::DataType)
	data_type, nothing
end
function parse_arg_info(p::Pair)
	parse_arg_info_pair(p.first, p.second)
end
function parse_arg_info(value)
	Any, value
end
function parse_arg_info_pair(a::DataType, b::DataType)
	if a isa b
		b, a
	elseif b isa a
		a, b
	else
		throw(ErrorException("DataTypeError: $(a) and $(b) are not inclusive."))
	end
end
parse_arg_info_pair(a, b::DataType) = b, convert_data_type(a, b)
parse_arg_info_pair(a::DataType, b) = a, convert_data_type(b, a)


function convert_data_type(value, data_type::DataType)
	if isa(value, data_type)
		return value
	else
		try
			convert(data_type, value)
		catch
			throw(ErrorException("DataTypeError: cannot convert $value to $data_type."))
		end
	end
end

## convert inputs/outputs to Dict{String} in `run(p, inputs, outputs)`
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
	for p in v
		res[str(p.first)] = p.second
	end
	res
end
to_xxput_dict(d::Dict{String}) = d
to_xxput_dict(any) = throw(ErrorException("DataTypeError: cannot run Program: cannot convert to Dict{String}: inputs, outputs, or returned value of infer_outputs(inputs)"))


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


## redirecting
Base.close(::Nothing) = nothing
Base.SimpleLogger(::Nothing) = nothing
Base.with_logger(f::Function, ::Nothing) = f()

Base.redirect_stdout(f::Function, ::Nothing) = f()
Base.redirect_stderr(f::Function, ::Nothing) = f()
Base.redirect_stdout(f::Function, ::Nothing) = f()

handle_open(::Nothing, mode) = nothing
handle_open(io::IO, mode) = io # do not change and do not close when exit
handle_open(file::AbstractString, mode) = open(file::AbstractString, mode)

"""
	redirect_to_files(f::Function, outfile; mode="a+")
	redirect_to_files(f::Function, outfile, errfile; mode="a+")
	redirect_to_files(f::Function, outfile, errfile, logfile; mode="a+")

Redirect outputs of function `f` to file(s).

- `xxxfile`: File path (`AbstractString`), `nothing` or `::IO`. `nothing` means no redirect. Files can be the same.
- `mode`: same as `open(..., mode)`.
"""
function redirect_to_files(f::Function, outfile, errfile, logfile; mode="a+")
	out = handle_open(outfile, mode)
	err = errfile == outfile ? out : handle_open(errfile, mode)
	log = logfile == outfile ? out : logfile == errfile ? err : handle_open(logfile, mode)
	res = redirect_stdout(out) do
		redirect_stderr(err) do
			logger = SimpleLogger(log)
			with_logger(logger) do
				f()
			end
		end
	end
	outfile isa IO || close(out)
	errfile isa IO || close(err)
	logfile isa IO || close(log)
	res
end



function redirect_to_files(f::Function, outfile::T, errfile::T; mode="a+") where T <: Union{Nothing, AbstractString}
	out = handle_open(outfile, mode)
	err = errfile == outfile ? out : handle_open(errfile, mode)
	res = redirect_stdout(out) do
		redirect_stderr(err) do
			logger = SimpleLogger(err)
			with_logger(logger) do
				f()
			end
		end
	end
	outfile isa IO || close(out)
	errfile isa IO || close(err)
	res
end

function redirect_to_files(f::Function, redirectfile::T; mode="a+") where T <: Union{Nothing, AbstractString}
	out = handle_open(redirectfile, mode)
	res = redirect_stdout(out) do
		redirect_stderr(out) do
			logger = SimpleLogger(out)
			with_logger(logger) do
				f()
			end
		end
	end
	redirectfile isa IO || close(out)
	res
end
