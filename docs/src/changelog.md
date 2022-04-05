# Change log

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
