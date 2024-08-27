# Change Log

## v0.11.0

- Feat/**Breaking**: new method `auto_change_directory(b::Bool)`. It is necessary because changing directory is not thread-safe in Julia. It was set to `false` in v0.11.0. To make your code compatible with previous version, you can add `Pipelines.auto_change_directory(true)` at the beginning of your code, or use full paths through out your code (recommended).

## v0.10.6

- Fix: `run(::Program)`: wrap `pwd()` in try-catch block in case the dir no longer exists. It happens because workding dir is not thread safe in Julia. If other program delete the directory, it will fail. 

## v0.10.5

- Compat: Julia v1.10: `hasmethod(f, t)` in Julia v1.10 changed its behavior. It affects Program function check. Use `length(methods(f, t))` instead.

## v0.10.4

- Compat: Julia v1.9: init(): change the if statement of printing Base.stdxxx was changed when initiating Pipelines.jl: now detects WindowsRawSocket.

## v0.10.3

- Compat: Julia v1.9: init(): change the if statement of printing Base.stdxxx was changed when initiating Pipelines.jl: now detects RawFD.

## v0.10.2

- Compat: Julia v1.9: init(): change the if statement of printing Base.stdxxx was changed when initiating Pipelines.jl: now detects `<fd>`
- Change `using FilePathsBase` to `import FilePathsBase:AbstractPath` to avoid variable conflict.
## v0.10.1

- Fix: the flaw of recursion of `keyword_interpolation`.

## v0.10.0

- Feature: remove reexport of FilePathsBase package to suppress warning 'Both Dates and FilePathsBase exports canonicalize'.

## v0.9.10

- Fix: `run`: wrap `cd(dir_backup)` in try-catch block in case the dir no longer exists. happens because workding dir is not thread safe in Julia. If other program delete the directory, it will fail.

## v0.9.9

- Fix: redirection.jl: wrap all `redirect_stdXXX` in the try-catch block.

## v0.9.8

- Feature: flush ios after task finished (`redirect_to_files`).

## v0.9.6 - v0.9.7

- Fix: not showing warning when Base.stdout or Base.stderr is system file redirection (IOStream, fd) when initiating Pipelines.jl

## v0.9.5

- Docs: add new page: tips and trobleshoots.

## v0.9.4

- Docs: better documentation on run id file: [`Pipelines.create_run_id_file`](@ref).
- Feature: [`Pipelines.CMD_FILE_SPLITER`](@ref).

## v0.9.3

- Feature: run id file also try to record all possible files from Base.AbstractCmd inputs/outputs.

## v0.9.2

- Feature: If skipping program, not showing `@info` of start.
- Feature: `@logmsg` does not show module, file, line etc.

## v0.9.1

- Fix: pretty print of Pipelines types.

## v0.9.0

- Feature: run id file contain file information to guess whether files are updated and better decide rerun.
- Fix: quote_function: func(a; b = b): if b is one of inputs or outputs, the first b does not change.

## v0.8.6

- Feature: If any file (not dir) paths of the inputs are newer than the run id file, the program will run again to update the outputs, even if the run id file exists. (#4)

## v0.8.5

- Feature: Enchance compatibility with JobSchedulers v0.6.12: `Program` has a new field called `arg_forward` that is used to forward user-defined inputs/outputs to specific keyword arguments of `JobSchedulers.Job(::Program, ...)`, including `name::String`, `user::String`, `ncpu::Int`, `mem::Int`.

## v0.8.4

- Fix: If defining a `Program` using `Expr` in other packages, the precompilation will fail due to the Expr will be `@eval` to a function in `Pipelines`. Julia does not allow define new functions in a compiled module. To avoid this, we have to pass the current module to `CmdProgram` or `JuliaProgram`. In this fix, the two methods accepts a new keyword argument `mod`, which allows users to pass `mod = @__MODULE__` manually.

## v0.8.2-3

- Feature: New `Arg` data type for storing inputs and outputs information in `Program`. `arg_inputs::Vector{Arg}` and `arg_outputs::Vector{Arg}` are new fields of `Program`.

- Feature: If a `Arg` name of inputs/outputs is a Symbol, run id will not generate using this Arg, which is useful for args do not affect the results, such as nthread, ncpu. (#8)

- Fix other issues.

- Docs: better documentation for v0.8 features.

## v0.8.1

- Optimize: `dictreplace!(ex::AbstractString, s::Symbol, v::Expr)` is now several times faster.

- Docs: better documentation for v0.8 features.

## v0.8.0

- Feature: `quote_function` allows user to use 'elements of inputs and outputs as variables' when defining Program. To use the feature, users can pass `quote ... end` to Program's arguments that requires `Function` before, such as `main`, `validate_inputs`. From Pipelines v0.8, all `Expr` provided will be converted to `Function` automatically. See details in `quote_expr`.

## v0.7.8

- Feature: `check_reserved_xxputs` does not need `xxput_types::Vector{Type}`.

## v0.7.7

- Feature: now `run(p; kwargs...)` works like `prog_run(p::Program; kwargs...)`. If original `inputs, outputs` are given, and same keys in `kwargs...` are found, the latter will override the former. If program args are conflict with other arguments, an error will throw.

- Feature: `infer_outputs(prog; INPUT1 = 5)` now supports `run`-like kwargs.

- Feature: Before creating new Program, `check_reserved_xxputs` and `check_function_methods`.

## v0.7.6

- Fix: replace `@run` in v0.7.5 with a new function `prog_run(p::Program; kwargs...)`. The original `@run` only works for global variables, otherwise you need to use `@eval` and complicated `$` to pass variables to AST.
- Fix replace `@vars` in v0.7.5 with a non-exported function `parse_program_args(p::Program; args...)`. It returns `(inputs::Dict{String}, outputs::Dict{String}, kwargs::Tuple)`.

## v0.7.5

- Feature (removed in v0.7.6): Simplify: new macro `@run program key_value_args... run_args...`: Run `program` without creating `inputs::Dict` and `outputs::Dict`. The inputs and outputs are provided in the form of `key = value`, rather than `Dict`.

- Feature (removed in v0.7.6): Simplify: new macro `@vars program::Program key_value_args...`: Return run-able `(inputs::Dict, outputs::Dict)` for `program` using `key_value_args` in the form of `key = value`.

## v0.7.4

- Fix: `check_dependency(p::Program)` return `::Bool` now.

- Feature: `check_dependency(p::CmdDependency; exit_when_fail::Bool = p.exit_when_fail)`: new argument `exit_when_fail` to override `p.exit_when_fail`.

- Feature: new function to check all `CmdDependency` and `Program` under `m::Module`: `check_dependency(m::Module = @__MODULE__; exit_when_fail = true, verbose = true)`

- Feature: new function to display dependency status: `status_dependency(m::Module = @__MODULE__; exit_when_fail = false, verbose = true)`.

## v0.7.3

- Fix a method overwrite warning when defining `Base.redirect_stdout(f::Function, ::Nothing) = f()`: redirect_stdout and redirect_stderr are the same (::Base.RedirectStdStream) at least from julia v1.7, so defining redirect_stdout means redirect_stderr is also defined. If diff exists in previous julia versions, check first.

## v0.7.2

- Optimize: do not show the error stack traces twice (`@error` and `rethrow()`).

## v0.7.1

- Optimize: the returned value of the main function of `JuliaProgram` does not required to be a `Dict`. If it is not a `Dict` or the returned Dict fails to pass keyword check, a warn message will be displayed and the inferred `outputs` will be used.

- Fix: MethodError: no method matching sort(::Dict{String, Any}) in `generate_run_uuid(...)`: import `OrderedCollections.jl`.

- Optimize: do not show the error messages twice when `stderr` is not redirected.

## v0.7.0

- Fix: `generate_run_uuid`: more stable way to generate run UUID using inputs and outputs.

**Caution**: The UUIDs generated using old versions will be out of date, so if you run the same Program with the same parameters, the program will generate a different UUID. This is why it is considered as a breaking version.

## v0.6.1

- Fix: `to_xxput_dict(d::Dict)`: `v` not defined: change `v` to `d`.

- Fix: `xxputs_completion_and_check(p, inputs, outputs)`: Now, if `outputs` are not empty and `p.infer_outputs` is defined, the function will first run `p.infer_outputs` and then merge the result and `outputs` (user-input keys are kept).

- Optimize: function `infer_outputs(p::Program, inputs)` now works like `xxputs_completion_and_check` but only return `outputs::Dict{String}`. It does not affect user-defined function `p.infer_outputs`. In addition, `inputs` is no longer strict to `Dict{String}` because the it will convert to `Dict{String}` first.

## v0.6.0

- Feature: allow retry for failed program: `run(p::Program, ...; retry = 1)`.

- Fix: the run id file generated after a successful run is under the directory specified in `run(p::Program, ...; dir = "directory")`, if `p.id_file` is not an absolute path (default).

- Feature: `verbose` in `run(p::Program, ...; verbose = :min)`: If `true` or `:all`, print all info and error messages. If `:min`, print minimum info and error messages. If `false` or `:none`, print error messages only.

## v0.5.2

- Fix: Export `restore_stderr`.

## v0.5.1

- Fix: stdxxx_origin set to nothing when reloading Pipelines from another module.

## v0.5.0

- Optimize: stack traces of failed results.

- Fix: original stdout and stderr can be recovered by using `restore_stdout()` and `restore_stderr()`. Remove `set_default_stdout()` and `set_default_stderr()` because they cause confusion sometimes.

- Optimize: default output of program info.

## v0.4.6

- Fix: results of `isok(::AbstractString)` should be expected.

## v0.4.5

- Optimize: allow check dependencies when test argument is empty.

## 0.4.4

- Fix: change DataType to Type throughout the code. It allows successful parsing of Union Type, such as `"READ1" => Union{String, Vector{String}}`.

## v0.4.3

- Fix file redirecting exception when redirecting to a closed stream. Not solved: Redirecting in Julia are not thread safe, so unexpected redirection might be happen if you are running programs in different `Tasks` or multi-thread mode.

## v0.4.2

- Better error stack trace after capturing.

- Update file redirecting.

## v0.4.1

- Fix file redirecting.

## v0.4.0

- Feature: `run(p::Program, ...)` supports running at a specified directory (`dir`). Run ID files will also create at that directory.

- Feature: `run(p::Program, ...)` supports redirecting `stdout`, `stderr` and `stdlog` (Julia log output, such as `@info`, `@warn`, `@error`).

- Fix: `CmdProgram` supports commands such as `pipeline(cmd, stdout=stderr)`, which redirect stdout to stderr now because `prepare_cmd(h::Base.TTY, inputs, outputs)` method is added.

## v0.3.2

- Fix: `CmdDependency`: do not check when `test_args` are empty.

## v0.3.1

- Feature: `to_cmd` and `to_str`: support argument `::Regex` or `::Any`.

## v0.3.0

- Feature: Building Program: Support type assertion and default arguments of `inputs` and `outputs`, such as `"arg" => 5`, `"arg" => Int`, `"arg" => 5 => Int`, `"arg" => Int => 5`.

- Feature: `Program` and `run(::Program)` no longer require `inputs` and `outputs` to be `Vector` or `Dict`. They can be both `Vector` or `Dict`, or even an element of `Vector` or `Dict`, as long as they can be converted. Eg:

  ```julia
  p = CmdProgram(
      cmd_dependencies = [julia],
      id_file = "id_file",
      inputs = [
          "input",
          "input2" => Int,
          "optional_arg" => 5,
          "optional_arg2" => 0.5 => Number
      ],
      outputs = "output" => "<input>.output"
      ,
      cmd = `echo input input2 optional_arg optional_arg2 output`
  )

  inputs = Dict(
      "input" => `in1`,
      "input2" => 2
  )

  outputs = [
      "output" => "out"
  ]

  run(p, inputs, outputs,
      skip_when_done = false,
      verbose = true,
      touch_run_id_file = false
  )
  ```

- Feature: Pretty print of `Program`. Eg:

  ```julia
  julia> p
  CmdProgram:
    name             → Command Program
    id_file          → id_file
    info_before      → auto
    info_after       → auto
    cmd_dependencies → CmdDependency[`/usr/software/julia-1.4.2/bin/julia -Cnative -J/usr/software/julia-1.4.2/lib/julia/sys.so -O3 -g1`]
    inputs           → "input"         :: Any    (required)
                       "input2"        :: Int64  (required)
                       "optional_arg"  :: Any    (default: 5)
                       "optional_arg2" :: Number (default: 0.5)
    validate_inputs  → do_nothing
    prerequisites    → do_nothing
    cmd              → `echo input input2 optional_arg optional_arg2 output`
    infer_outputs    → do_nothing
    outputs          → "output" :: Any (default: <input>.output)
    validate_outputs → do_nothing
    wrap_up          → do_nothing
  ```

## v0.2.2

- Support JobSchedulers.jl.

## v0.2.1

- Fix examples in docs.

## v0.2.0

- CmdDependency: Better interpolation in `Cmd`.

  ```julia
  dep::CmdDepencendy

  # old version
  `$(dep.exec) --args`
  # or
  `$(exec(dep)) --args`

  # now
  `$dep --args`
  ```

- New `JuliaProgram` for pure Julia implementation.

- `Program` is the Abstract type containing `CmdProgram` and `JuliaProgram` substypes.
