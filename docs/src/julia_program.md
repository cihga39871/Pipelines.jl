# Julia Program

Please read the Command Program section in Manual page first.

Pipelines are built with multiple `Program`s. `Program` is the abstract type contains `CmdProgram` and `JuliaProgram`.

## Differences between CmdProgram and JuliaProgram

`JuliaProgram` and `CmdProgram` are generally the same and remain most compatibility, except for the differences:

| Diff                 | `cp :: CmdProgram`                                           | `jp :: JuliaProgram`                                         |
| :------------------- | :----------------------------------------------------------- | :----------------------------------------------------------- |
| Unique Field         | `cp` has a command template: `cp.cmd::AbstractCmd`           | `jp` has a main Julia function: `jp.main::Function`          |
| Method to Run        | Replace keywords in `cp.cmd` with values in user defined `inputs` and `outputs::Dict{String}` | Simply evoke `jp.main(inputs, outputs::Dict{String})`        |
| **Returned Outputs** | `outputs ` **provided in** `run(...)` **does not change**    | `outputs` **is overwritten by the returned value of** `jp.main` |
| Dry Run              | Return `(replaced_cmd::AbstractCmd, run_id_file::String)`    | Return `(fake_outputs::Dict{String}, run_id_file::String)`   |

## Structure

`JuliaProgram` can be built with this method:

```julia
JuliaProgram <: Program

JuliaProgram(;
	name::String               = "Unnamed Command Program",
	id_file::String            = "",
	info_before::String        = "auto",
	info_after::String         = "auto",
	cmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),
	inputs                     = Vector{String}(),
	validate_inputs::Function  = do_nothing,  # positional arguments: inputs::Dict{String}
	prerequisites::Function    = do_nothing,  # positional arguments: inputs, outputs::Dict{String}
	main::Function             = do_nothing,  # positional arguments: inputs, outputs::Dict{String}
	infer_outputs::Function    = do_nothing,  # positional arguments: inputs::Dict{String}
	outputs                    = Vector{String}(),
	validate_outputs::Function = do_nothing  # positional arguments: outputs::Dict{String},
	wrap_up::Function          = do_nothing  # positional arguments: inputs, outputs::Dict{String}
) -> JuliaProgram
```

## Run

To run a `JuliaProgram`, the methods are the same as `CmdProgram`:

```julia
run(
	p::JuliaProgram;
	inputs=Dict{String}(),
	outputs=Dict{String}(),
	skip_when_done::Bool=true,
	check_dependencies::Bool=true,
	stdout=nothing,
	stderr=nothing,
	append::Bool=false,
	verbose::Bool=true,
	touch_run_id_file::Bool=true,
	dry_run::Bool=false
) -> (success::Bool, outputs::Dict{String})

run(p::JuliaProgram, inputs, outputs; kwargs...)

run(p::JuliaProgram, inputs; kwargs...)
)  # only usable when `p.infer_outputs` is defined, or default outputs are set in `p`.
```

!!! note "Compatibility with JobSchedulers.jl"

    Pipelines.jl is fully compatible with [JobSchedulers.jl](https://github.com/cihga39871/JobSchedulers.jl) which is a Julia-based job scheduler and workload manager inspired by Slurm and PBS.
    
    `run(::Program, ...)` can be replaced by `Job(::Program, ...)`. The latter creates a `Job`, and you can submit the job to queue by using `submit!(::Job)`. See example below.

## Example
```julia
using Pipelines

p = JuliaProgram(
	id_file = "id_file",
	inputs = ["a", "b"],
	outputs = ["out"],
	main = (inputs, outputs) -> begin
		a = inputs["a"]
		b = inputs["b"]
		println("inputs are ", a, " and ", b)
		println("You can also use info in outputs: ", outputs["out"])

		println("The returned value will be assigned to a new outputs")
		return Dict{String,Any}("out" => b^2)
	end
)

inputs = Dict(
	"a" => `in1`,
	"b" => 2
)

outputs = Dict(
	"out" => "will_be_replaced"
)

success, outputs = run(p, inputs, outputs;
	touch_run_id_file = false
) # outputs will be refreshed
```

### Compatibility with JobSchedulers.jl

```julia
using JobSchedulers

scheduler_start()  # start job scheduler

job = Job(p, inputs, outputs;
	touch_run_id_file = false
)  # create a Job object; same arguments as `run`

submit!(job)  # submit job to queue

result(job)  # get the returned result

```
