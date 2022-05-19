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

## Quote for Program Args
```@docs
quote_expr
quote_function
```

## Run Programs

```@docs
Base.run
```

## Command Dependency
```@docs
CmdDependency
check_dependency
status_dependency
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

## Internal
```@docs
RESERVED_KEY_SET
parse_default
to_xxput_dict
try_function
StackTraceVector
keyword_interpolation
xxputs_completion_and_check
parse_program_args
```
