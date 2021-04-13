# Julia Program

Please read the Command Program section in Manual page first.

Pipelines are built with multiple `Program`s. `Program` is the abstract type contains `CmdProgram` and `JuliaProgram`.

## Differences between CmdProgram

`JuliaProgram` and `CmdProgram` are generally the same and remain most compatibility, except for the differences:

- A `cp::CmdProgram` contains a command template (defined in `cp.cmd`), while a `jp::JuliaProgram` contains a main Julia function (defined in `jp.main`).
- The main Julia function takes `inputs::Dict{String}` and `outputs::Dict{String}` as positional arguments. That means when we call `run(jp::JuliaProgram, ...)`, `jp.main(inputs::Dict{String}, outputs::Dict{String})` will be evaluated.
- The returned value of the main function will be assigned to new `outputs`, ie. it is evaluated like this `outputs = jp.main(inputs, outputs)`. Please ensure the returned value is `Dict{String}` with proper keys.
- When using dry run (`run(jp::JuliaProgram, dry_run=true)`), return `(fake_outputs::Dict{String}, run_id_file::String)`.

## Structure

`JuliaProgram` can be built with this method:

```julia
JuliaProgram(;
	name::String               = "Julia Program",
	id_file::String            = "",
	info_before::String        = "auto",
	info_after::String         = "auto",
	cmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),
	inputs::Vector{String}     = Vector{String}(),
	validate_inputs::Function  = do_nothing,  # positional arguments: inputs::Dict{String}
	prerequisites::Function    = do_nothing,  # positional arguments: inputs, outputs::Dict{String}
	main::Function             = do_nothing,  # positional arguments: inputs, outputs::Dict{String},
	infer_outputs::Function    = do_nothing,  # positional arguments: inputs::Dict{String}
	outputs::Vector{String}    = Vector{String}(),
	validate_outputs::Function = do_nothing,  # positional arguments: outputs::Dict{String}
	wrap_up::Function          = do_nothing  # positional arguments: inputs, outputs::Dict{String}
)
```

To run a `JuliaProgram`, the methods are the same as `CmdProgram`:

```julia
run(
	p::JuliaProgram;
	inputs::Dict{String}=Dict{String, Cmd}(),
	outputs::Dict{String}=Dict{String, Cmd}(),
	skip_when_done::Bool=true,
	check_dependencies::Bool=true,
	stdout=nothing,
	stderr=nothing,
	append::Bool=false,
	verbose::Bool=true,
	touch_run_id_file::Bool=true,
	dry_run::Bool=false
) -> (success::Bool, outputs::Dict{String})

run(
	p::JuliaProgram,
	inputs::Dict{String},
	outputs::Dict{String};
	kwargs...
)

run(
	p::JuliaProgram,
	inputs::Dict{String},
	kwargs...
)  # only usable when `p.infer_outputs` is defined.
```

## Example
```julia
p = JuliaProgram(
	cmd_dependencies = [julia],
	id_file = "id_file",
	inputs = ["a", "b"],
	outputs = ["c"],
	main = (inputs, outputs) -> begin
		println("inputs are ", inputs["a"], " and ", inputs["b"])
		println("You can also use info in outputs:", outputs["c"])
        println("The returned value will be assigned to a new outputs")
        return Dict{String,Any}("c" => b^2)
	end
)

inputs = Dict(
	"a" => `in1`,
	"b" => 2
)

outputs = Dict(
	"c" => "out"
)

success, outputs = run(p, inputs, outputs;
	touch_run_id_file = false
) # outputs will be refreshed
```
