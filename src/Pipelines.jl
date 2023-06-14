module Pipelines

using Reexport
using Printf
@reexport using Dates
using UUIDs
using Logging
using OrderedCollections  # sorting dict in generate_run_uuid
import FilePathsBase:AbstractPath  # isinputnewer in Program.jl

include("quote_function.jl")
export quote_function, quote_expr

include("utils.jl")
export do_nothing, isok,
str, to_str, to_cmd,
replaceext, removeext

include("redirection.jl")
export restore_stdout, restore_stderr,
redirect_to_files

include("CmdDependency.jl")
export CmdDependency,
check_dependency,
check_dependency_file,
check_dependency_dir

include("Arg.jl")
export Arg

include("Program.jl")
export Program, prog_run,
infer_outputs,
status_dependency

include("run_id_file.jl")


include("CmdProgram.jl")
export CmdProgram,
prepare_cmd

include("JuliaProgram.jl")
export JuliaProgram

include("pretty_print.jl")

# binding documentations
@doc (@doc prog_run) run

function __init__()
    global stdout_origin
    global stderr_origin

    if isnothing(stdout_origin)
        if Base.stdout isa Base.TTY
            nothing
        elseif occursin(r"<fd .*>|RawFD\(\d+\)", string(Base.stdout))
            nothing
        else
            # Not Terminal (TTY), nor linux file redirection (fd)
            @warn "Base.stdout was changed when initiating Pipelines.jl. `restore_stdout()` can only restore to the current one." Base.stdout
        end
        stdout_origin = Base.stdout
    end
    if isnothing(stderr_origin)
        if Base.stderr isa Base.TTY
           nothing 
        elseif occursin(r"<fd .*>|RawFD\(\d+\)", string(Base.stderr))
            nothing
        else
            # Not Terminal (TTY), nor linux file redirection (fd)
            @warn "Base.stderr was changed when initiating Pipelines.jl. `restore_stderr()` can only restore to the current one." Base.stderr
        end
        stderr_origin = Base.stderr
    end
end

end # module
