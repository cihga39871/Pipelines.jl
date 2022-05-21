module Pipelines

using Dates
using UUIDs
using Logging
using OrderedCollections  # sorting dict in generate_run_uuid

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
        Base.stdout isa Base.TTY || (@warn "Base.stdout was changed when initiating Pipelines.jl. `restore_stdout()` can only restore to the current one." Base.stdout)
        stdout_origin = Base.stdout
    end
    if isnothing(stderr_origin)
        Base.stderr isa Base.TTY || (@warn "Base.stderr was changed when initiating Pipelines.jl. `restore_stderr()` can only restore to the current one." Base.stderr)
        stderr_origin = Base.stderr
    end
end

end # module
