# API

## Command Program

```@docs
CmdProgram()
run(::CmdProgram)
infer_outputs(::CmdProgram, ::Dict{String})
```

## Command Dependency

```@docs
CmdDependency()
check_dependency(::CmdDependency)
check_dependency_dir(path::AbstractString; exit_when_false=true)
check_dependency_file(path::AbstractString; exit_when_false=true)
```

## Utils

```@docs
to_str(::Cmd)
to_cmd(::Cmd)
split(::Cmd)
replaceext(::String, ::AbstractString)
removeext(::String)
```
