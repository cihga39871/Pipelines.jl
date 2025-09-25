module Pipelines

using Reexport
using Printf
@reexport using Dates
using UUIDs
using Logging
using OrderedCollections  # sorting dict in generate_run_uuid
import FilePathsBase:AbstractPath  # isinputnewer in Program.jl

@reexport using ScopedStreams

const redirect_to_files = ScopedStreams.redirect_stream
export redirect_to_files

include("quote_function.jl")
export quote_function, quote_expr

include("utils.jl")
export do_nothing, isok,
str, to_str, to_cmd,
replaceext, removeext

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

ScopedStreams.@gen_scoped_stream_methods true

end # module
