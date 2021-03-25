# API

## Command Dependency

```@docs
CmdDependency()
check_dependency(::CmdDependency)
check_dependency_dir(path::AbstractString; exit_when_false=true)
check_dependency_file(path::AbstractString; exit_when_false=true)
```

## Command Program

```@docs
CmdProgram()
run(::CmdProgram)
infer_outputs(::CmdProgram, ::Dict{String})
```