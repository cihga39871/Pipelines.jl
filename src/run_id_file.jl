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

Line format of run id file (tab delimited):

- Column 1: `i` or `o` stands for inputs or outputs.
- Column 2: file modified time in Float64.
- Column 3: file size.
- Column 3: key name.
- Column 4: file path.
"""
function create_run_id_file(run_id_file::AbstractString, inputs::Dict, outputs::Dict)
    tmp_file = tmpname()
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
    mv(string(tmp_file), run_id_file; force=true)
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

Check whether re-run the program. Return `true` means need re-run.
"""
function need_rerun(p::Program, run_id_file::AbstractString, inputs::Dict, outputs::Dict)
    isfile(run_id_file) || (return true)

    isok(getfield(p, :validate_outputs)(outputs)) || (return true)

    any_file_differ(run_id_file::AbstractString, inputs::Dict, outputs::Dict) && (return true)

    false
end

"""
    RUN_ID_LINE_SKIP_EXTENSION = $RUN_ID_LINE_SKIP_EXTENSION

If a file with a extension listed, infomation of this file will not write to `run_id_file`.
"""
const RUN_ID_LINE_SKIP_EXTENSION = [".so", ".dylib", ".dll"]

function write_run_id_line(io, arg_name::AbstractString, file::AbstractString, first_char::String; skip_ext::Vector = RUN_ID_LINE_SKIP_EXTENSION)
    ext = splitext(file)[2]
    ext in skip_ext && (return false)  # skip dynamic library

    modified_time = mtime(file)
    fs = filesize(file)
    println(io, "$first_char\t$(modified_time)\t$fs\t$arg_name\t$file")
    true
end
write_run_id_line(io, arg_name::AbstractPath, file::AbstractString, first_char::String; skip_ext::Vector = RUN_ID_LINE_SKIP_EXTENSION) = write_run_id_line(io, string(arg_name), file, first_char; skip_ext = skip_ext)


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
                splited = split(f, r"[,;:]")
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
            splited = split(arg, r"[,;:]")
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
