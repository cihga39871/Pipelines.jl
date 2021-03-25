# Command Program

Please read the Quick Start section in Home page first.

Pipelines are built with multiple `CmdProgram`s. A `CmdProgram` contains a command template, lists of dependencies/inputs/outputs, and predefined methods to prepare run-time environment, validate dependencies/inputs/outputs, and so on.

## Structure

`CmdProgram` can be build with the method:

```julia
CmdProgram(;
	name::String               = "Unnamed Command Program",
	id_file::String            = "",
	info_before::String        = "auto",
	info_after::String         = "auto",
	cmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),
	inputs::Vector{String}     = Vector{String}(),
	validate_inputs::Function  = do_nothing,  # positional arguments: inputs::Dict{String}
	prerequisites::Function    = do_nothing,  # positional arguments: inputs, outputs::Dict{String}
	cmd::Base.AbstractCmd      = ``,
	infer_outputs::Function    = do_nothing,  # positional arguments: inputs::Dict{String}
	outputs::Vector{String}    = Vector{String}(),
	validate_outputs::Function = do_nothing,  # positional arguments: outputs::Dict{String},
	wrap_up::Function          = do_nothing   # positional arguments: inputs, outputs::Dict{String}
) -> CmdProgram
```

In this way, all preparation and post-evaluation can be wrapped in a single `CmdProgram`. It is easy to maintain and use.

To run a `CmdProgram`, use one of the following methods:

```julia
run(
	p::CmdProgram;
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
	p::CmdProgram,
	inputs::Dict{String},
	outputs::Dict{String};
	kwargs...
)

run(
	p::CmdProgram,
	inputs::Dict{String},
	kwargs...
)  # only usable when `p.infer_outputs` is defined.
```

The explanation of arguments is in the next section.

## Workflow

1. **Check keywords consistency:** Inputs/outputs keywords should be consistent in both `p::CmdProgram` and `run(p; inputs, outputs)`.

   > For example, if inputs and outputs in `p::CmdProgram` is defined this way
   >
   > ```julia
   > p = CmdProgram(..., inputs = ["I", "J"], outputs = ["K"])
   > ```
   >
   > The inputs and outputs in `run(...)` has to be
   >
   > ```julia
   > run(p;
   >     inputs = Dict(
   >         "I" => something,
   >         "J" => something
   >     ),
   >     outputs = Dict("K" => something)
   > )
   > ```

2. Print info about starting program. 

   > The content can set by `p.info_before::String`.
   >
   > Disable: `run(..., verbose=false)`

3. Check whether the program ran successfully before. If so, return `(true, outputs::Dict{String})` without running it.

   !!! note "How does the program know it ran before?"
       a. `run(..., skip_when_down=true)` enables the functionality.
       
   
       b. Run id file is generated. The prefix of run id file is set by `p.id_file`. After given inputs and outputs, a unique ID will be append to the prefix. You can use `run(..., touch_run_id_file=false)` to skip creating the run id file.
       
       c. `p.validate_outputs(outputs)` run successfully without returning `false`.
   
4. Check command dependencies (`CmdDependency`). 
   
> Disable: `run(..., check_dependencies=false)`
   >
   > Read **Command Dependency** portion for details.
   
5. Remove the run id file if exists.

6. Validate inputs by invoking `p.validate_inputs(inputs)`.

7. Preparing the main command.

   > If you specify `run(...; stdout=something, stderr=something, append::Bool)`, the command (`cmd`) will be wrapped with `pipeline(cmd; stdout=something, stderr=something, append::Bool)`. If `cmd` has its own file redirection, the outer wrapper may not work as you expect.

8. Meet prerequisites by invoking `p.prerequisites(inputs, outputs)`.

   > It is the last function before running the main command.  For example, you can create a directory if  the main command cannot create itself.

9. Run the main command.

10. Validate outputs by invoking `p.validate_outputs(outputs)`.

11. Run the wrap up function by invoking `p.wrap_up(inputs, outputs)`

    > It is the last function to do after-command jobs. For example, you can delete intermediate files if necessary.

12. Create run id file if `run(..., touch_run_id_file=true)`. Read Step 3 for details.

13. Print info about finishing program. 

    > The content can set by `p.info_after::String`.
    >
    > Disable: `run(..., verbose=false)`

14. Return `(success::Bool, outputs{String})`