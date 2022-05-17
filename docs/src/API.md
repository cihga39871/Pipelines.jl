# API

## Program
```@docs
Program
infer_outputs
parse_program_args
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

### New way (since v0.7.7)
```@docs
prog_run
```

### Old way
```@docs
Base.run
```

## Command Dependency
```@docs
CmdDependency
check_dependency
check_dependency_dir(path::AbstractString; exit_when_false=true)
check_dependency_file(path::AbstractString; exit_when_false=true)
```

## Common Methods
```@docs
status_dependency
```

## Utils
```@docs
replaceext(::String, ::AbstractString)
removeext(::String)
split(::Cmd)
to_str(::Cmd)
to_cmd(::Cmd)
do_nothing
isok
```

## Redirection
```@docs
redirect_to_files
restore_stdout
restore_stderr
```
