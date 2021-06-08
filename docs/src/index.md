# Pipelines.jl

```@meta
CurrentModule = Pipelines
```

*A lightweight Julia package for computational pipelines.*

Building reusable pipelines and workflows is easier than you have ever thought.

## Package Features

- Easy to build both simple and complex tasks.

- Supports external command lines and pure Julia functions.

- Supports **resuming** interrupted tasks, **skipping** finished tasks.

- Supports dependency check.

- Supports inputs, outputs validation, and so on.

- Supports program queuing and workload management with [JobSchedulers.jl](https://github.com/cihga39871/JobSchedulers.jl)

## Installation

Pipelines.jl can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run

```julia
pkg> add Pipelines
```

To use the package, type

```julia
using Pipelines
```

## Quick Start

Pipelines are built with multiple `Program`s. `Program` is the abstract type of `CmdProgram` and `JuliaProgram`.

A `CmdProgram` contains a command template and name lists of inputs/outputs. The names of inputs/outputs will be replaced by real values when executing the program.

Let's set up a simple `CmdProgram` to print values using `echo`:

```julia
using Pipelines

echo = CmdProgram(
    inputs = [
        "REQUIRED",               # no default value; any data type.
        "TYPED" => String,        # no default value; String type only.
        "OPTIONAL" => 4,          # default value is 5; any data type.
        "FULL" => String => "abc" # default value is abc; String type only.
    ],
    cmd = `echo REQUIRED TYPED OPTIONAL FULL`   
)
```

Running the program is just like running other `Cmd`,  but here we need to specify inputs by using `Dict{String => value}` (`Vector{String => value}` is also supported.)

```julia
inputs = Dict(
    "REQUIRED" => "Pipelines.jl",
    "TYPED" => "is",
    "FULL" => "everyone!"
)
run(echo, inputs)
```

!!! note "Program will not run twice by default!"
    If you run a program with the same inputs again, the program will just return the same result, display a warning message without running the command twice.

    ```julia
    run(echo, inputs)
    ```

    This is because the program will generate a file (run id file) in the current directory indicating the program has been run. Several methods can be used to re-run a program:

    ```julia
    # Method 1: stop checking finished program
    run(echo, inputs; skip_when_done = false)

    # Method 2: delete the run_id_file before running again:
    cmd, run_id_file = run(echo, inputs; dry_run = true) # Dry-run returns the command and run id file without running it.
    rm(run_id_file)  # remove the run_id_file

    # Method 3: Do not generate run_id_file when first running.
    run(echo, inputs; touch_run_id_file=false)
    ```

### Program with Outputs

Unlike the first example, many programs write files as outputs. Pipelines.jl has an elegant way to handle it.

The following program prints values simultaneously, sort them, and save to a file.

```julia
prog = CmdProgram(
    inputs = ["INPUT1", "INPUT2", "INPUT3"],
    outputs = "OUTPUT_FILE",
    cmd = pipeline(`echo INPUT1 INPUT2` & `echo INPUT3`, `sort`, "OUTPUT_FILE")
)

inputs = Dict(
    "INPUT1" => "Hello,",
    "INPUT2" => `Pipeline.jl`,
    "INPUT3" => 39871
)
outputs = "OUTPUT_FILE" => "out.txt" # save output to file

run(prog, inputs, outputs) # will return (success::Bool, outputs)

run(`cat out.txt`) # print the content of out.txt
# 39871
# Hello, Pipeline.jl
```

### Default values

Default values and data types can be set for keywords of `inputs` and `outputs` in this way:

```julia
echo = CmdProgram(
    inputs = [
        "REQUIRED",                     # no default value; any data type.
        "TYPED" => String,              # no default value; String type only.
        "OPTIONAL" => 4,                # default value is 5; any data type.
        "FULL1" => String => "abc"      # default value is abc; String type only.
        "FULL2" => "abc" => String      # default value is abc; String type only.
        "INTERPOLATED" => "<FULL1>.xyz" # default value is value of FULL1 * ".xyz".
    ],
    cmd = `echo REQUIRED TYPED OPTIONAL FULL`   
)
```

#### Interpolation of default values

If the default value is a `String`, it can be interpolated by using `<keyword>`, such as `"<FULL1>.xyz"` in the example.

#### Generate outputs using Function

> This step is prior to adding default values of outputs, and string interpolation using `<>`.

We also provide a method (`infer_outputs::Function`) in `CmdProgram` to generate complex `outputs::Dict{String}` from `inputs::Dict{String}`. The type used in the function is restricted to `Dict{String}`

```julia
using Dates

prog = CmdProgram(
    inputs = [
        "INPUT1" => Int,
        "INPUT2" => Int => 3
    ],
    outputs = "OUTPUT_FILE",
    cmd = pipeline(`echo INPUT1 INPUT2`, `sort`, "OUTPUT_FILE"),
    infer_outputs = inputs -> Dict(
    	"OUTPUT_FILE" => string(now(), "__", inputs["INPUT1"], ".txt")
    )
)
success, outputs = run(prog, "INPUT1" => 5)
```

We can also generate default outputs without running the program:

```julia
outputs = infer_outputs(prog, inputs)
```

### Julia Program

Pipelines also defined `JuliaProgram` type for pure Julia functions. It is like `CmdProgram` and remain most compatibility. More details are in the Julia Program, Manual Page.

### Compatibility with JobSchedulers.jl

Pipelines.jl is fully compatible with [JobSchedulers.jl](https://github.com/cihga39871/JobSchedulers.jl) which is a Julia-based job scheduler and workload manager inspired by Slurm and PBS.

`run(::Program, ...)` can be replaced by `Job(::Program, ...)`. The latter creates a `Job`, and you can submit the job to queue by using `submit!(::Job)`.

## Future Development

- Support running competitive tasks with **locks**.

## Change log

v0.3.0

- Building Program: Support type assertion and default arguments of `inputs` and `outputs`, such as `"arg" => 5`, `"arg" => Int`, `"arg" => 5 => Int`, `"arg" => Int => 5`.

- `Program` and `run(::Program)` no longer require `inputs` and `outputs` to be `Vector` or `Dict`. They can be both `Vector` or `Dict`, or even an element of `Vector` or `Dict`, as long as they can be converted. Eg:

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

- Pretty print of `Program`. Eg:

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

v0.2.2

- Support JobSchedulers.jl.

v0.2.1

- Fix examples in docs.

v0.2.0

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
