module Pipelines

using Logging, LoggingExtras, Dates
using UUIDs

include("utils.jl")
export do_nothing, isok,
str, to_str, to_cmd,
replaceext, removeext


include("CmdDependency.jl")
export CmdDependency,
check_dependency,
check_dependency_file,
check_dependency_dir,
exec

include("CmdProgram.jl")
export CmdProgram,
infer_outputs



end # module
