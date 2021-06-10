module Pipelines

using Dates
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

include("Program.jl")
export Program,
infer_outputs

include("CmdProgram.jl")
export CmdProgram,
prepare_cmd

include("JuliaProgram.jl")
export JuliaProgram

include("pretty_print.jl")

end # module
