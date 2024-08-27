function generate_run_uuid(p::Program, inputs::Dict{String}, outputs::Dict{String})
    out_uuid = UUID4
    in_names = sort!([arg.name for arg in p.arg_inputs])
    out_names = sort!([arg.name for arg in p.arg_outputs])
    for name in in_names
        out_uuid = uuid5(out_uuid, string(name, ":", inputs[name]))
    end
    for name in out_names
        out_uuid = uuid5(out_uuid, string(name, ":", outputs[name]))
    end
    return out_uuid
end

function run_id_file_path(p::Program, dir::AbstractString, run_id)
    joinpath(dir, getfield(p, :id_file) * "." * string(run_id))
end

# """
#     isinputnewer(inputs::Dict, run_id_file::AbstractString) -> Bool

# Check whether any existing file (not dir) paths of `AbstractString` or `AbstractPath` in `inputs` are newer than `run_id_file`.
# """
# function isinputnewer(inputs::Dict, run_id_file::AbstractString)
#     input_time = 0.0
#     for (k,v) in inputs
#         if v isa AbstractString || v isa AbstractPath
#             if isfile(v)
#                 modified_time = mtime(v)
#                 if modified_time > input_time
#                     input_time = modified_time
#                 end
#             end
#         end
#     end

#     input_time > mtime(run_id_file)
# end

mutable struct RunIDLine
    isinput::Bool
    mtime::Float64
    filesize::Int
    key::String
    value::String
end
function RunIDLine(line::AbstractString)
    vals = split(chomp(line), '\t')
    if length(vals) != 5
        @error "Cannot parse this line of run id file: $line"
        return nothing
    end
    try
        RunIDLine(
            vals[1] == "i",
            parse(Float64, vals[2]),
            parse(Int, vals[3]),
            vals[4],
            vals[5]
        )
    catch e
        @error "Cannot parse this line of run id file: $line"
        return nothing
    end
end

"""
    create_run_id_file(run_id_file::AbstractString, inputs::Dict, outputs::Dict)

Create run id file.

## What is run id file

The run id file stores information of arguments and files related to a successful run of `Program`.

By comparing the run id file, we can determine whether we need to re-run a finished program.

## File Name

The name of a run id file is `<dir>/<program_id_prefix>.<argument_UUID>`.

- `dir`: working directory to run the program, which can be defined in `run(::Program; dir = "")`

- `program_id_prefix`: the prefix of run ID file, which can be defined in `CmdProgram(id_file = "")` and `JuliaProgram(id_file = "")`

- `argument_UUID`: a unique ID generated from string representations of inputs and outputs arguments using the internal function `generate_run_uuid`.

In this way, the name of a run id file will not change if running a program in the same directory with same inputs and outputs.

However, this is not enough for determine whether a job needs re-run. Consider this situation:

> (1) Run prog with arg = 1, output `"out.txt"` and `"run_id_file_with_arg1"`  
>
> (2) Run prog with arg = 2, output `"out.txt"` and `"run_id_file_with_arg2"`  
>
> (3) Run prog with arg = 1 again, no re-run because `"out.txt"` and `"run_id_file_with_arg1"` all exist!  

To solve the issue, we need to store the states of inputs and outputs arguments.

Here, we guess file names from inputs and outputs. 

- If an argument is `AbstractString` or `AbstractPath`, and `isfile()` returns true, we store the file information. (We ignore directories because their contents are easy to change.)

- If an argument is `Base.AbstractCmd`, we decompose the command into pieces, and check whether each piece is a file path. The rules of file guessing are complicated and mentioned in [`Pipelines.cmd_to_run_id_lines`](@ref) and [`Pipelines.CMD_FILE_SPLITER`](@ref).

- If a file name is found, and its extension is not one of [`Pipelines.RUN_ID_LINE_SKIP_EXTENSION`](@ref), it will write to run id file.

!!! note "Control file name guessing"
    Pipeline developers usually know what file extension should be ignored, and whether they have an argument joining two files with a splitter. In this way, we can use IN-PLACE methods to change [`Pipelines.RUN_ID_LINE_SKIP_EXTENSION`](@ref) and [`Pipelines.CMD_FILE_SPLITER`](@ref). IN-PLACE methods are usually functions ending with `!`, such as `empty!`, `push!`, `deleteat!`

## Contents of run id file

- Tab delimited, no header.
- Column 1: `i` or `o` stands for inputs or outputs.
- Column 2: unix timestamp of when the file was last modified in Float64.
- Column 3: the size (in bytes) of the file.
- Column 3: key name of inputs or outputs. It may have duplication.
- Column 4: file path. It may have duplication.

## Limitation

We cannot store states of all arguments. If we have a pure JuliaProgram without reading and writing files, we cannot guarantee the state of the arguments.

A work-around is to intentionally create a file with a fixed name, and the file name is defined in Program's outputs.

## See also
[`Pipelines.cmd_to_run_id_lines`](@ref), [`Pipelines.RUN_ID_LINE_SKIP_EXTENSION`](@ref), [`Pipelines.CMD_FILE_SPLITER`](@ref)
"""
function create_run_id_file(run_id_file::AbstractString, inputs::Dict, outputs::Dict)
    tmp_file = tempname()
    open(tmp_file, "w+") do io
        for (k,v) in inputs
            if v isa AbstractString || v isa AbstractPath
                if isfile(v)
                    write_run_id_line(io, k, v, "i")
                end
            elseif v isa Base.AbstractCmd
                cmd_to_run_id_lines(io, k, v, "i")
            end
        end

        for (k,v) in outputs
            if v isa AbstractString || v isa AbstractPath
                if isfile(v)
                    write_run_id_line(io, k, v, "o")
                end
            elseif v isa Base.AbstractCmd
                cmd_to_run_id_lines(io, k, v, "o")
            end
        end
    end
    try
        mv(tmp_file, run_id_file; force=true)
    catch
        @warn "Cannot create run id file at $run_id_file. Skip creating."
        isfile(tmp_file) && rm(tmp_file)
    end
end

function parse_run_id_file(run_id_file::AbstractString)
    res = Dict{String, RunIDLine}()
    for line in eachline(run_id_file)
        run_id_line = RunIDLine(line)
        isnothing(run_id_line) && continue
        res[run_id_line.value] = run_id_line
    end
    res
end

"""
    any_file_differ(run_id_file::AbstractString, inputs::Dict, outputs::Dict)

Check whether any existing file (not dir) path of `AbstractString` or `AbstractPath` differ from records in `run_id_file`.
"""
function any_file_differ(run_id_file::AbstractString, inputs::Dict, outputs::Dict)
    file_info = parse_run_id_file(run_id_file)
    any_file_differ(file_info, inputs) && (return true)
    any_file_differ(file_info, outputs) && (return true)
    false
end

function any_file_differ(file_info::Dict{String, RunIDLine}, xxputs::Dict)
    for (k,v) in xxputs
        if v isa AbstractString || v isa AbstractPath
            if isfile(v)
                v = string(v)
                now_mtime = mtime(v)
                old = get(file_info, v, nothing)
                if old === nothing
                    return true
                else
                    (now_mtime == old.mtime) || (return true)
                    (filesize(v) == old.filesize) || (return true)
                end
            else # file missing, if it is input, ignore. outputs = differ
                old = get(file_info, v, nothing)
                old === nothing && continue
                old.isinput || (return true)
            end
        elseif v isa Base.AbstractCmd
            # scan file_info with the same arg_name
            # if old has any change, re run
            for old in values(file_info)
                if old.key == k
                    f = old.value
                    isfile(f) || (return true)
                    (mtime(f) == old.mtime) || (return true)
                    (filesize(f) == old.filesize) || (return true)
                end
            end
        end
    end
    false
end

"""
    need_rerun(p::Program, run_id_file::AbstractString, inputs::Dict, outputs::Dict) -> Bool

Check whether re-run the program `p`. Return `true` means it need re-run.

## Decision details

1. Is `run_id_file` a file? If not, re-run.

2. Run `p.validate_outputs(outputs)`. If fail, re-run.

3. Comparing status of files in `run_id_file` using [`Pipelines.any_file_differ`](@ref). If yes, re-run.
"""
function need_rerun(p::Program, run_id_file::AbstractString, inputs::Dict, outputs::Dict)
    isfile(run_id_file) || (return true)

    isok(getfield(p, :validate_outputs)(outputs)) || (return true)

    any_file_differ(run_id_file::AbstractString, inputs::Dict, outputs::Dict) && (return true)

    false
end

"""
    RUN_ID_LINE_SKIP_EXTENSION = $RUN_ID_LINE_SKIP_EXTENSION

If a file with an extension listed, `run_id_file` skip storing information of this file. It means whether to re-run a program, the state of the file will be ignored.

See also: [`Pipelines.create_run_id_file`](@ref), [`Pipelines.cmd_to_run_id_lines`](@ref), [`Pipelines.CMD_FILE_SPLITER`](@ref)

"""
const RUN_ID_LINE_SKIP_EXTENSION = String[".so", ".dylib", ".dll"]

"""
    CMD_FILE_SPLITER = $CMD_FILE_SPLITER

It is aimed to guess whether an argument of a command contain multiple file names joined using file splitters. 

See also: [`Pipelines.create_run_id_file`](@ref), [`Pipelines.cmd_to_run_id_lines`](@ref), [`Pipelines.RUN_ID_LINE_SKIP_EXTENSION`](@ref)
"""
const CMD_FILE_SPLITER = Char[',', ';', ':']

function write_run_id_line(io, arg_name::AbstractString, file::AbstractString, first_char::String; skip_ext::Vector = RUN_ID_LINE_SKIP_EXTENSION)
    ext = splitext(file)[2]
    ext in skip_ext && (return false)  # skip dynamic library

    modified_time = mtime(file)
    fs = filesize(file)
    println(io, "$first_char\t$(modified_time)\t$fs\t$arg_name\t$file")
    true
end
write_run_id_line(io, arg_name::AbstractPath, file::AbstractString, first_char::String; skip_ext::Vector = RUN_ID_LINE_SKIP_EXTENSION) = write_run_id_line(io, string(arg_name), file, first_char; skip_ext = skip_ext)

"""
    cmd_to_run_id_lines(io::IO, arg_name::AbstractString, cmd::Base.AbstractCmd, first_char::String)

- `io`: IO of run id file.
- `arg_name`: name of the inputs/outputs argument.
- `cmd`: Subtypes of Base.AbstractCmd.
- `first_char`: `"i"` or `"o"`, stands for inputs or outputs.

## Rules to guess file names from command:

- The first argument is ignored because usually it is a script.

- Numbers are ignored.

- If an argument starts with `-`, matching `r"^-[A-Za-z0-9\\-\\_]+(/.+)"` and `r"^-[A-Za-z0-9\\-\\_]+=(.+)"` only. If matched, go to the next rule.

- Check whether an arg is a file. If not, try to use `Pipelines.CMD_FILE_SPLITER` to split the argument, and check each part. If found a file, go to the next rule.

- If a file name is found, and its extension is not one of `Pipelines.RUN_ID_LINE_SKIP_EXTENSION`, it will write to run id file.
"""
function cmd_to_run_id_lines(io::IO, arg_name::AbstractString, cmd::Cmd, first_char::String)
    for (i, arg) in enumerate(cmd.exec)
        i == 1 && continue  # the first argument is software, skip
        length(arg) == 0 && continue  # empty argument
        tryparse(Float64, arg) isa Float64 && continue  # not file

        if arg[1] == '-' && length(arg) > 2
            # like -qJ/mnt/scratch_2t/usr/software/julia-1.8.1/lib/julia/sys.so
            m = match(r"^-[A-Za-z0-9\-\_]+(/.+)", arg)

            if isnothing(m)
                m = match(r"^-[A-Za-z0-9\-\_]+=(.+)", arg)
                if isnothing(m)
                    continue
                end
            end

            f = m.captures[1]
            if isfile(f)
                write_run_id_line(io, arg_name, f, first_char)
            else
                # check whether it is file1,file2 or file1;file2
                splited = split(f, CMD_FILE_SPLITER)
                length(splited) == 1 && continue
                for s in splited
                    if isfile(s)
                        write_run_id_line(io, arg_name, s, first_char)
                    end
                end
            end
            continue
        end

        if isfile(arg)
            write_run_id_line(io, arg_name, arg, first_char)
        else
            # check whether it is file1,file2 or file1;file2
            splited = split(arg, CMD_FILE_SPLITER)
            length(splited) == 1 && continue
            for f in splited
                if isfile(f)
                    write_run_id_line(io, arg_name, f, first_char)
                end
            end
        end
    end
end

function cmd_to_run_id_lines(io::IO, arg_name::AbstractString, c::Base.CmdRedirect, first_char::String)
    cmd_to_run_id_lines(io, arg_name, c.cmd, first_char)
    cmd_to_run_id_lines(io, arg_name, c.handle, first_char)
end
function cmd_to_run_id_lines(io::IO, arg_name::AbstractString, c::T, first_char::String) where T <: Union{Base.OrCmds, Base.ErrOrCmds, Base.AndCmds}
    cmd_to_run_id_lines(io, arg_name, c.a, first_char)
    cmd_to_run_id_lines(io, arg_name, c.b, first_char)
end
function cmd_to_run_id_lines(io::IO, arg_name::AbstractString, h::Base.TTY, first_char::String)   # h: handle property of CmdRedirect
    nothing
end
function cmd_to_run_id_lines(io::IO, arg_name::AbstractString, h::Base.FileRedirect, first_char::String)   # h: handle property of CmdRedirect
    if isfile(h.filename)
        write_run_id_line(io, arg_name, h.filename, first_char)
    end
end
