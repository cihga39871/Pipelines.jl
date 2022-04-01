# API

## Program
```@docs
Program
```

## Command Program
```@docs
CmdProgram()
prepare_cmd
```

## Julia Program
```@docs
JuliaProgram()
```

## Run Programs
```@docs
Base.run
```

## Common Methods
```@docs
infer_outputs(::Program, ::Dict{String})
check_dependency(::Program)
```

## Command Dependency
```@docs
CmdDependency()
check_dependency(::CmdDependency)
check_dependency_dir(path::AbstractString; exit_when_false=true)
check_dependency_file(path::AbstractString; exit_when_false=true)
exec
```

## Utils
```@docs
to_str(::Cmd)
to_cmd(::Cmd)
split(::Cmd)
replaceext(::String, ::AbstractString)
removeext(::String)
```

## Redirection
```@docs
redirect_to_files
restore_stdout
restore_stderr
```
