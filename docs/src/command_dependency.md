# Command Dependency

Command Dependency (`CmdDependency`) is the external dependent program that can be used in Command Program (`CmdProgram`).

## Create

```@docs
CmdDependency()
```

## Check

`CmdDependency` can be checked manually by using

```julia
check_dependency(p::CmdDependency)
```

Usually, it is not necessary because program builder will add necessary `CmdDependency` to `CmdProgram`, and dependency check will be done when running `CmdProgram` by default.

## Use with `CmdProgram`

When building `CmdProgram`, we can specify our command dependencies by using

```julia
dep1 = CmdDependency(...)
dep2 = CmdDependency(...)

prog = CmdProgram(
	cmd_dependencies = [dep1, dep2],
    ...
)
```

To call the dependency in command, we can use `exec(dep)` or `dep.exec`:

```julia
prog = CmdProgram(
	cmd_dependencies = [dep1, dep2],
    cmd = `$(exec(dep1)) --args` & `$(dep2.exec) --args`
)
```

### An Example

We use `samtools` and `bowtie2` as dependencies when building a mapping program.

```julia
SAMTOOLS = CmdDependency(
	exec = `samtools`,
	test_args = `--version`,
	validate_success = true,
	validate_stdout = x -> occursin(r"^samtools \d+", x)
)

BOWTIE2 = CmdDependency(
	exec = `bowtie2`,
	test_args = `--version`,
	validate_success = true,
	validate_stdout = x -> occursin(r"bowtie2-align-s version \d+", x)
)

program_bowtie2 = CmdProgram(
    name = "Bowtie2 Mapping",
    id_file = ".bowtie2",
    cmd_dependencies = [BOWTIE2, SAMTOOLS],
    
    inputs = ["FASTQ", "REF"],
    validate_inputs = inputs -> begin
        check_dependency_file(inputs["FASTQ"]) && 
        check_dependency_file(inputs["REF"])
    end,
    
    outputs = ["BAM"],
	infer_outputs = inputs -> begin
        Dict("BAM" => str(inputs["FASTQ"]) * ".bam")
    end,
    
    cmd = pipeline(`$(exec(BOWTIE2)) -x REF -q FASTQ`, `$(exec(SAMTOOLS)) sort -O bam -o BAM`),
    
    wrap_up = (inputs, outputs) -> begin
        run(`$(exec(SAMTOOLS)) index $(outputs["BAM"])`)
    end
)
```

To run the program, we can decide whether to check dependencies with `run(..., check_dependencies=true)`:

```julia
inputs = Dict("FASTQ" => "a.fastq", "REF" => "ref.fasta")
run(program_bowtie2, inputs; check_dependencies = true)
```

