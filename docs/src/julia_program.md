# Julia Program

Please read the Command Program section in Manual page first.

Pipelines are built with multiple `Program`s. `Program` is the abstract type contains `CmdProgram` and `JuliaProgram`.

## Differences between CmdProgram and JuliaProgram

`JuliaProgram` and `CmdProgram` are generally the same and remain most compatibility, except for the differences:

| Diff                 | `cp :: CmdProgram`                                           | `jp :: JuliaProgram`                                         |
| :------------------- | :----------------------------------------------------------- | :----------------------------------------------------------- |
| Unique Field         | `cp` has a command template: `cp.cmd::AbstractCmd`           | `jp` has a main Julia code: `jp.main` that can be defined using `quote .. end` . |
| Method to Run        | Replace keywords in `cp.cmd` with values in user defined `inputs` and `outputs}` | If `main` is a `Expr`, the variables mentioned in inputs and outputs will be replaced to the `inputs["VAR"]` format, and then a Julia function returning `outputs::Dict{String}` is created and called. If `main` isa `Function`, it just invokes `jp.main(inputs::Dict{String}, outputs::Dict)`. |
| **Returned Outputs** | `outputs ` **provided in** `run(...)` **does not change**    | `outputs` **will only be overwritten by the returned value of** `jp.main` when the returned value is a `Dict` and passes `p.validate_outputs`. |
| Dry Run              | Return `(replaced_cmd::AbstractCmd, run_id_file::String)`    | Return `(fake_outputs::Dict{String}, run_id_file::String)`   |

## Structure

`JuliaProgram` can be built with this method:

```julia
JuliaProgram <: Program

JuliaProgram(;
    name::String                            = "Julia Program",
    id_file::String                         = "",
    info_before::String                     = "auto",
    info_after::String                      = "auto",
    cmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),
    inputs                                  = Vector{String}(),
    validate_inputs::Expr                   = do_nothing,  # vars of inputs
    infer_outputs::Expr                     = do_nothing,  # vars of inputs
    prerequisites::Expr                     = do_nothing,  # vars of inputs and outputs
    main::Expr                              = do_nothing,  # vars of inputs and outputs
    outputs                                 = Vector{String}(),
    validate_outputs::Expr                  = do_nothing,  # vars of outputs
    wrap_up::Expr                           = do_nothing   # vars of inputs and outputs
) -> JuliaProgram
```

## Run

To run a `JuliaProgram`, the methods are the same as `CmdProgram`:

```julia
success, outputs = run(
	p::Program;
	program_kwargs...,
	dir::AbstractString="",
	check_dependencies::Bool=true,
	skip_when_done::Bool=true,
	touch_run_id_file::Bool=true,
	verbose=true,
	retry::Int=0,
	dry_run::Bool=false,
	stdout=nothing,
	stderr=nothing,
	stdlog=nothing,
	append::Bool=false
) -> (success::Bool, outputs::Dict{String})
```

- `program_kwargs...` include elements in `p.inputs` and `p.outputs`
- Other keyword arguments are related to run. Details can be found at [`run`](@ref).



!!! warning "Thread safety"
    Redirecting and directory change in Julia are not thread safe, so unexpected redirection and directory change might be happen if you are running programs in different `Tasks` or multi-thread mode.



!!! compat "Compatibility with JobSchedulers.jl"

    Pipelines.jl is fully compatible with [JobSchedulers.jl](https://github.com/cihga39871/JobSchedulers.jl) which is a Julia-based job scheduler and workload manager inspired by Slurm and PBS.

    `run(::Program, ...)` can be replaced by `Job(::Program, ...)`. The latter creates a `Job`, and you can submit the job to queue by using `submit!(::Job)`. See example below.

## Example
```julia
using Pipelines

p = JuliaProgram(
	id_file = "id_file",
	inputs = ["a",
	          "b" => Int],
	outputs = "c" => "<a>.<b>",
	main = quote
		println("inputs are ", a, " and ", b)
		println("You can also use info in outputs: ", outputs["c"])
		println("The returned value will be assigned to a new outputs")
		println("It is ok to use inputs and outputs directly:")
		@show inputs
		@show outputs
		c = b^2
	end)

# running the program using `run`: keyword arguments include keys of inputs and outputs
success, new_out = run(p; a = `in1`, b = 2, c = "out", touch_run_id_file = false)

# an old way to `run` program: need to create inputs and outputs first.
inputs = Dict("a" => `in1`, "b" => 2)
outputs = "c" => "out"
success, new_out = run(p, inputs, outputs; touch_run_id_file = false)

# for CmdProgram, outputs are inferred before running the main command, however,
# for JuliaProgram, outputs will change to the returned value of main function, if the returned value is a Dict and pass `p.validate_outputs`
@assert new_out != outputs
```

### Compatibility with JobSchedulers.jl

```julia
using JobSchedulers

scheduler_start()  # start job scheduler

job = Job(p, a=`in1`, b=2, out="any", touch_run_id_file=false)

submit!(job)  # submit job to queue

result(job)  # get the returned result

```
