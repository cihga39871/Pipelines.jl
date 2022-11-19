# API

## Program
```@docs
Program
infer_outputs
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

## Quote for Program
```@docs
quote_expr
quote_function
```

## Arg
```@docs
Arg
Pipelines.RESERVED_KEY_SET
Pipelines.FORWARD_KEY_SET
```

## Run Program
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

## Internal
```@docs
Pipelines.parse_arg
Pipelines.to_xxput_dict
Pipelines.try_function
Pipelines.StackTraceVector
Pipelines.keyword_interpolation
Pipelines.xxputs_completion_and_check
Pipelines.parse_program_args
Pipelines.RUN_ID_LINE_SKIP_EXTENSION
Pipelines.create_run_id_file
Pipelines.CMD_FILE_SPLITER
Pipelines.cmd_to_run_id_lines
Pipelines.need_rerun
Pipelines.any_file_differ
```
