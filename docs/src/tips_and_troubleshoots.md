# Tips & Troubleshoots

This page summarizes tips and troubleshoots when building workflow using Pipelines.jl

## Named tuple for common arguments

A workflow usually has many `Program`s, and all `Programs` have some common run arguments. We can use a named tuple to store and use the common arguments easily.

```julia
using Pipelines

prog_A = CmdProgram(...)
prog_B = JuliaProgram(...)

common_run_args = (check_dependencies = false, skip_when_done = true, touch_run_id_file = true, verbose = :min, retry = 1)

run(prog_A; prog_A_args..., common_run_args...)
run(prog_B; prog_B_args..., common_run_args...)
```

## Run different programs in parallel

When building a computational workflow, we may find different programs use different CPU and memory. Some can run simultaneously, but some have to run sequentially. To efficiently use computational resources, we highly recommend to use [JobSchedulers.jl](https://github.com/cihga39871/JobSchedulers.jl). It is stable, useful and powerful package for task queuing and workload management, and are fully compatible with Pipelines.jl

`run(prog; prog_args..., run_args...)` can be replaced by the JobScheduler way:

```julia
job = Job(prog; prog_args..., run_args..., job_args...)
submit!(job)
```

### An example

A workflow comprises the following steps:

- Run `prog_A` with 2 threads and 4GB RAM.
- Run `prog_B` with 8 threads.
- After `prog_A` finished, run `prog_C` (2 threads).
- After `prog_B` and `prog_C` finished, run `prog_D` (12 threads)

```julia
using JobSchedulers, Pipelines

scheduler_start() # start the job scheduler

prog_A = CmdProgram(...)
prog_B = JuliaProgram(...)
prog_C = CmdProgram(...)
prog_D = JuliaProgram(...)

job_A = Job(prog_A, A_args..., ncpu = 2, mem = 4GB)
submit!(job_A)

job_B = Job(prog_B, B_args..., ncpu = 8)
submit!(job_B)

job_C = Job(prog_C, C_args..., ncpu = 2, dependency = DONE => job_A)
submit!(job_C)

job_D = Job(prog_D, D_args..., ncpu = 12, dependency = [DONE => job_B, DONE => job_C])
submit!(job_D)
```

### Argument forward

Usually a parallel program has an input argument of # cpu, and `Job()` has a `ncpu` argument. Because `ncpu` is reserved ([`Pipelines.RESERVED_KEY_SET`](@ref)) and cannot be used as Program inputs, it is awkward to write code in this way:

```julia
prog = CmdProgram(inputs = [:NCPU => 1], ...)  # default num CPU is 1
job = Job(prog, NCPU = 8, ncpu = 8)            # use 8 CPUs in the job
```

To make life easier, we can define `arg_forward` in `Program` to forward the input `:NCPU` to Job argument `ncpu`:

```julia
prog = CmdProgram(inputs = [:NCPU => 1], arg_forward = :NCPU => :ncpu, ...)
job = Job(prog, NCPU = 8)
# Job:
#   ...
#   ncpu → 8
#   ...
```

!!! note "arg_forward in Program"
    `arg_forward`: forward args from inputs and outputs to specific keywords in `JobSchedulers.Job()`, only supporting `Pipelines.FORWARD_KEY_SET`: `[:ncpu, :mem, :user, :name]`. Elements (or vectors containing elements) in the following format: `"arg_of_inputs_or_outputs" => :key_in_FORWARD_KEY_SET`.

## Name of inputs and outputs: String or Symbol?

Normally the name should be a String: `CmdProgram(inputs = ["IN"], outputs = ["OUT"])`. However, if an argument does not affect results (such as number of threads), it is called an *independent* argument, and has to be a **Symbol**. Symbol arguments are **ignored** when generating unique run IDs to prevent **re-running** a program. Arguments of inputs and outputs will be converted to [`Arg`](@ref) objects.

See also: [`Arg`](@ref)


## UndefVarError in `quote ... end`

`quote` creates a piece of code without using the explicit `Expr` constructor. It creates an `Expr` object and follows the scoping rules of Julia `Expr`ession.

An `Expr`ression will be evaluated in the **global scope** of the **Module** defined in `Program(..., mod = Module)`, *and then converted to a function using [`Pipelines.quote_function`](@ref)*.

Thus, directly using local variables (including functions) in `Expr` will lead to `UndefVarError`. To use a local variable, you need to follow the rules:

1. A local variable (include function) should be referenced using `$` in `expr`ession.

2. A local `::Symbol` variable (`sym`) should be referenced using `$(QuoteNode(sym))` in `expr`ession.

Example:

```julia
inputs = ["A", "B"]
g_var = 3
g_sym = :globalsymbol

function gen_expr()
    l_var = 5
    l_func() = @info("Use local function")
    l_sym = :abc
    expr = quote
        @show inputs
        @show g_var
        @show g_sym
        @show $(QuoteNode(l_sym))
        @show $l_var + 2
        $l_func()
        A + B
    end
end

expr = gen_expr()
func = Pipelines.quote_function(expr, inputs; mod = @__MODULE__)

in_dict = Dict("A" => 5, "B" => 50) # func takes Dict{String} as argument
func(in_dict)
```

See also: [`quote_expr`](@ref)

## Fail to Precompile a module containing a `Program`

It is because you pass `quote ... end` to a Program, but forget to add `mod = @__MODULE__` to it.

Fix is simple:

```julia
prog = CmdProgram(...)                    # fail precompilation
prog = CmdProgram(..., mod = @__MODULE__) # pass precompilation
```

See also: [`quote_expr`](@ref)

## I changed `Program` inputs, but `Program` skipped and outputs not updated

The `Program` probably does not read or write files, and you may adjust inputs and outputs back and forth.

If we have a pure Program without reading and writing files, we cannot guarantee the state of the arguments.

A work-around is to intentionally create a file with a fixed name, and the file name is defined in Program's outputs.

See [`Pipelines.create_run_id_file`](@ref) to learn how Pipelines.jl decides whether a program needs re-run, and its limitation.
