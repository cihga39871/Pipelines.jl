

# logging with dates
const date_format = "yyyy-mm-dd HH:MM:SS"
timestamp_logger(logger) = TransformerLogger(logger) do log
    merge(log, (; message = "$(Dates.format(now(), date_format)) $(log.message)"))
end
ConsoleLogger(stdout, Logging.Debug) |> timestamp_logger |> global_logger

# UUID generation
const UUID4 = uuid4(UUIDs.MersenneTwister(39871))

# funcs
do_nothing() = nothing
do_nothing(x) = nothing
do_nothing(x, y) = nothing

"""
    str(c::Cmd) = string(c)[2:end-1]

Convert `Cmd` to `String` and remove backticks.
"""
function str(c::Cmd)
    string(c)[2:end-1]
end
to_str = str

to_cmd(c::Cmd) = c
to_cmd(s::AbstractString) = Cmd([string(s)])
to_cmd(v::Vector{String}) = Cmd(v)
to_cmd(v::Vector{T}) where T = Cmd(string.(v))
to_cmd(::Nothing) = nothing
to_cmd(i::Number) = Cmd([string(i)])

ValidInputTypes = Union{Cmd, AbstractString, Vector, Number}

"""
    split(c::Cmd) = c.exec :: Vector{String}

Return splitted arguments of `Cmd`.
"""
function Base.split(c::Cmd)
    c.exec
end

isok(x::Nothing) = true
isok(x::Bool) = x
isok(x::AbstractString) = occursin(r"y(es)?|ok?|t(rue)?|^1$"i, x)
