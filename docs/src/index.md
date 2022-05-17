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
        "OPTIONAL" => 4,          # default value is 4; any data type.
        "FULL" => String => "abc" # default value is abc; String type only.
    ],
    cmd = `echo REQUIRED TYPED OPTIONAL FULL`   
)
```

Running the program is just like running other `Cmd`,  but here we need to specify inputs in keyward arguments. Caution: use `;` to split positional and keyward arguments, and do not use `,`.

```julia
prog_run(echo; REQUIRED = "Pipelines", TYPED = "are", FULL = "to build.", OPTIONAL = :easy)
```

!!! note "Program will not run twice by default!"
    If you run a program with the same inputs again, the program will just return the same result, display a warning message without running the command twice.

    ```julia
    prog_run(echo; REQUIRED = "Pipelines", TYPED = "are", FULL = "to build.", OPTIONAL = :easy)
    ```

    This is because the program will generate a file (run id file) in the current directory indicating the program has been run. Several methods can be used to re-run a program:

    ```julia
    # Method 1: stop checking finished program using skip_when_done = false
    prog_run(echo; input_args..., skip_when_done = false)

    # Method 2: delete the run_id_file before running again
    cmd, run_id_file = prog_run(echo; input_args..., dry_run = true) # Dry-run returns the command and run id file without running it.
    rm(run_id_file)  # remove the run_id_file

    # Method 3: Do not generate run_id_file after a successful run
    prog_run(echo; input_args..., touch_run_id_file=false)
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

# inputs and outputs can be mixed together:
prog_run(prog; INPUT1 = "Good", INPUT2 = "Morning", INPUT3 = "Human", OUTPUT_FILE = "morning.txt")

run(`cat morning.txt`) # print the content of out.txt
# Good Morning
# Human
```

### Default values

Default values and data types can be set for keywords of `inputs` and `outputs` in this way:

```julia
echo = CmdProgram(
    inputs = [
        "REQUIRED",                     # no default value; any data type.
        "TYPED" => String,              # no default value; String type only.
        "OPTIONAL" => 4,                # default value is 5; any data type.
        "FULL1" => String => "abc",     # default value is abc; String type only.
        "FULL2" => "abc" => String,     # default value is abc; String type only.
        "INTERPOLATED" => "<FULL1>.xyz" # default value is value of FULL1 * ".xyz".
    ],
    cmd = `echo REQUIRED TYPED OPTIONAL FULL`   
)
```

#### Interpolation of default values

If the default value is a `String`, it can be interpolated by using `<keyword>`, such as `"<FULL1>.xyz"` in the example.

#### Generate outputs using Function

> This step is prior to adding default values of outputs, and string interpolation using `<>`.

We also provide a parameter (`infer_outputs::Function`) in `CmdProgram` to generate complex `outputs::Dict{String}` from `inputs::Dict{String}`. The argument (`inputs`) and returned value of the function has to be a `Dict{String}`.

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
success, outputs = prog_run(prog; INPUT1 = 5)
```

We can also generate default outputs without running the program:

```julia
inputs = Dict("INPUT1" => 5)
outputs = infer_outputs(prog, inputs)
```

### Julia Program

Pipelines also defined `JuliaProgram` type for pure Julia functions. It is like `CmdProgram` and remain most compatibility. More details are in the Julia Program, Manual Page.

### Compatibility with JobSchedulers.jl

Pipelines.jl is fully compatible with [JobSchedulers.jl](https://github.com/cihga39871/JobSchedulers.jl) which is a Julia-based job scheduler and workload manager inspired by Slurm and PBS.

`prog_run(::Program, ...)` can be replaced by `Job(::Program, ...)`. The latter creates a `Job`, and you can submit the job to queue by using `submit!(::Job)`.

## Future Development

- Support running competitive tasks with **locks**.
