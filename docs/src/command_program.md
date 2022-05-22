# Command Program

Please read the Quick Start section in Home page first.

Pipelines are built with multiple `Program`s. `Program` is the abstract type contains `CmdProgram` and `JuliaProgram`.

A `CmdProgram` contains a command template, lists of dependencies/inputs/outputs, and predefined methods to prepare run-time environment, validate dependencies/inputs/outputs, and so on.


## An Example

We will go through an example to illustrate how to write a `CmdProgram`.

The example is a **robust** Bowtie2 mapping program. It shows every feature of `CmdProgram`, but in reality, we do not need that much validation.

### Define the Name
```julia
program_bowtie2 = CmdProgram(
    name = "Bowtie2 Mapping",# The name of CmdProgram
    id_file = ".bowtie2"     # When job completed, a file ".bowtie2.xxxxxx" will
                             # be created to indicate the job is finished to
                             # avoid re-run.
)
```

### Main Command (Required)

In Bash script, the main code of the example:

```bash
REF=/path/to/reference_genome
FASTQ=/path/to/input_fastq_file
BAM=/path/to/output_bam_file
NTHREADS=8

bowtie2 --threads $NTHREADS -x $REF -q $FASTQ | samtools sort -O bam -o $BAM
```

The equivalent version using CmdProgram:

```julia
program_bowtie2 = CmdProgram(
    ...,
    inputs = [
        "FASTQ" => String,
        "REF" => "human_genome_hg38.fa" => String,
        :NTHREADS => Int => 8
    ],
    outputs = "BAM" => String,
    cmd = pipeline(`bowtie2 --threads NTHREADS -x REF -q FASTQ`, `samtools sort -O bam -o BAM`)
)
```

Now, the code can be run by invoking `run(program_bowtie2; FASTQ = "x", REF = "y", NTHREADS = 8, BAM = "z")`, but to illustrate all features, we will add more things to make it robust and easy to use.

!!! tip "Name of arguments: String or Symbol"
    Normally the name should be a `String`. However, if the argument does not affect results (such as number of threads), it needs to be a `Symbol`. Symbol arguments are ignored when generating unique run IDs to prevent re-running a program. Arguments will be converted to [`Arg`](@ref) objects.

### Command Dependency (Robustness↑)

We use `samtools` and `bowtie2` (two bioinformatics programs) as `command dependencies`. They can be wrapped in [`CmdDependency`](@ref).

### Infer Outputs (Convenience↑)

We can set a default method to generate `outputs::Dict{String}` from inputs, which allows us run the program without specifying outputs. Elements in `inputs` can be directly used as variables.

```julia
program_bowtie2 = CmdProgram(
    ...,
    outputs = "BAM" => String,
    infer_outputs = quote
        Dict("BAM" => FASTQ * ".bam")
    end,
    ...
)
```

!!! tip "Quote usage in Program"
    Quote is a piece of code of the type `Expr`ession. See details in [`quote_expr`](@ref).

Or, the following does the same job:

```julia
program_bowtie2 = CmdProgram(
    ...,
    outputs = "BAM" => "<FASTQ>.bam" => String,
    ...
)
```

!!! tip "Useful functions to change file names"
    - [`replaceext`](@ref): replace file extension.
    - [`removeext`](@ref): remove file extension.
    - [`to_str`](@ref): converts most types to `String`.
    - [`to_cmd`](@ref): converts most types to `Cmd`.
    More details are in the API/Utils page.


### Validate Inputs (Robustness↑)

To make the code robust, we can check whether the inputs exists by using `validate_inputs`. Elements in `inputs` can be directly used as variables.

```julia
program_bowtie2 = CmdProgram(
    ...,
    inputs = [
        "FASTQ" => String,
        "REF" => "human_genome_hg38.fa" => String
    ],
    validate_inputs = quote
        check_dependency_file(FASTQ) && check_dependency_file(REF)
    end,
    ...
)
```

### Prerequisites (Robustness↑)

We also need to prepare something (`prerequisites`) before running the main command. For example, create the output directory if not exist. Elements in `inputs` and `outputs` can be directly used as variables.

```julia
program_bowtie2 = CmdProgram(
    ...,
    prerequisites = quote
        mkpath(dirname(BAM))
    end
)
```

### Validate Outputs (Robustness↑)

After running the main command, we can validate outputs by using `validate_outputs`. Elements in `outputs` can be directly used as variables.

```julia
program_bowtie2 = CmdProgram(
    ...,
    outputs = "BAM",
    validate_outputs = quote
        check_dependency_file(BAM)
    end,
    ...
)
```

### Wrap Up (Convenience↑)

After validating outputs, we may also do something to wrap up, such as removing temporary files. Here, we build an index for output BAM file. Elements in `inputs` and `outputs` can be directly used as variables.

```julia
program_bowtie2 = CmdProgram(
    ...,
    wrap_up = quote
        run(`samtools index $BAM`)  # dollar sign is necessary in quote, unlike Pipelines(;cmd = ...) cannot use dollar sign.
    end
)
```

### The Final Code

All in all, the final program is like this:

```julia
program_bowtie2 = CmdProgram(
    name = "Bowtie2 Mapping",
    id_file = ".bowtie2",

    inputs = [
        "FASTQ" => String,
        "REF" => "human_genome_hg38.fa" => String
    ],

    outputs = ["BAM" => String],
    infer_outputs = quote
        Dict("BAM" => FASTQ * ".bam")
    end,

    validate_inputs = quote
        check_dependency_file(FASTQ) && check_dependency_file(REF)
    end,

    prerequisites = quote
        mkpath(dirname(BAM))
    end,

    validate_outputs = quote
        check_dependency_file(BAM)
    end,

    cmd = pipeline(`bowtie2 -x REF -q FASTQ`, `samtools sort -O bam -o BAM`),  # do not use dollar sign here.

    wrap_up = quote
        run(`samtools index $BAM`)  # unlike cmd = ..., dollar sign is necessary in all quotes!
    end
)
```

## Structure

`CmdProgram` can be built with this method:

```julia
CmdProgram <: Program

CmdProgram(;
    name::String                            = "Command Program",
    id_file::String                         = "",
    info_before::String                     = "auto",
    info_after::String                      = "auto",
    cmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),
    inputs                                  = Vector{String}(),
    validate_inputs::Expr                   = do_nothing,  # vars of inputs
    infer_outputs::Expr                     = do_nothing,  # vars of inputs
    prerequisites::Expr                     = do_nothing,  # vars of inputs and outputs
    cmd::Base.AbstractCmd                   = ``,
    outputs                                 = Vector{String}(),
    validate_outputs::Expr                  = do_nothing,  # vars of outputs
    wrap_up::Expr                           = do_nothing   # vars of inputs and outputs
) -> CmdProgram
```

In this way, all preparation and post-evaluation can be wrapped in a single `CmdProgram`. It is easy to maintain and use.

## Run

To run a `Program`, use this method:

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


!!! compat
    The old methods < v0.8 still work, which put program's arguments in `inputs::Dict{String}` and `outputs::Dict{String}`:

    ```julia
    success, outputs = run(p::Program, inputs, outputs; run_kwargs...)
    # only usable when outputs have default values.
    success, outputs = run(p::Program, inputs; run_kwargs...)
    ```



!!! warning
    Redirecting and directory change in Julia are not thread safe, so unexpected redirection and directory change might be happen if you are running programs in different `Tasks` or multi-thread mode.



!!! compat "Compatibility with JobSchedulers.jl"

    Pipelines.jl is fully compatible with [JobSchedulers.jl](https://github.com/cihga39871/JobSchedulers.jl) which is a Julia-based job scheduler and workload manager inspired by Slurm and PBS.

    `run(::Program, ...)` can be replaced by `Job(::Program, ...)`. The latter creates a `Job`, and you can submit the job to queue by using `submit!(::Job)`.

The explanation of arguments is in the next section.

## Workflow

1. Go to the working directory. Establish redirection. (`dir`, `stdout`, `stderr`, `stdlog`, `append`).

2. **Check keywords consistency:** Inputs/outputs keywords should be consistent in both `p::CmdProgram` and `run(p; inputs..., outputs...)`.

   > For example, if inputs and outputs in `p::CmdProgram` is defined this way
   >
   > ```julia
   > p = CmdProgram(..., inputs = ["I", "J"], outputs = ["K"])
   > ```
   >
   > You have to provide all I, J and K:
   >
   > ```julia
   > run(p; I = something, J = something, K = something)
   > ```

3. Print info about starting program.

   > The content can set by `p.info_before::String`.
   >
   > Disable: `run(..., verbose=false)`

4. Check whether the program ran successfully before. If so, return `(true, outputs::Dict{String})` without running it.

   > **How does the program know it ran before?**
   >
   > a. `run(..., skip_when_done=true)` skip running the program if it has been done before.
   >
   > b. Run id file is generated. The prefix of run id file is set by `p.id_file`. After given inputs and outputs, a unique ID will be appended to the prefix. You can use `run(..., touch_run_id_file=false)` to skip creating the run id file.
   >
   > c. `p.validate_outputs(outputs)` run successfully without returning `false`.

5. Check command dependencies (`CmdDependency`).

   > Disable: `run(..., check_dependencies=false)`
   >
   > Read **Command Dependency** portion for details.

6. Remove the run id file if exists.

7. Validate inputs. (`p.validate_inputs`)

8. Preparing the main command.

   > If you specify `run(...; stdout=something, stderr=something, append::Bool)`, the command (`cmd`) will be wrapped with `pipeline(cmd; stdout=something, stderr=something, append::Bool)`. If `cmd` has its own file redirection, the outer wrapper may not work as you expect.

9. Meet prerequisites. (`p.prerequisites`)

   > It is the last code before running the main command.  For example, you can create a directory if  the main command cannot create itself.

10. Run the main command.

11. Validate outputs. (`p.validate_outputs`)

12. Run the wrap up code. (`p.wrap_up`)

    > It is the last code to do after-command jobs. For example, you can delete intermediate files if necessary.

13. Create run id file if `run(..., touch_run_id_file=true)`. Read Step 3 for details.

14. Print info about finishing program.

    > The content can set by `p.info_after::String`.
    >
    > Disable: `run(..., verbose=false)`
    >
    > Simple info: `run(..., verbose=min)`

15. Return `(success::Bool, outputs{String})`

!!! note "Dry Run"
    `run(..., dry_run=true)` will return `(mature_command::AbstractCmd, run_id_file::String)` instead.
