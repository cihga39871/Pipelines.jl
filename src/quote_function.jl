
"""
    quote
        do_some_thing()
    end :: Expr

`quote` creates a piece of code without using the explicit `Expr` constructor.

From Pipelines v0.8, you can use `quote ... end` to `validate_inputs`, `infer_outputs`, do `prerequisites`, do `main`, `validate_outputs`, and `wrap_up` a `Program`.

Elements in `inputs` or `outputs` can be directly used as variables for those arguments. See the table below.

| Argument                 | Elements as variables | Default returned value                             |
| :----------------------- | :-------------------- | :------------------------------------------------- |
| validate_inputs          | inputs                | the last expression                                |
| infer_outputs            | inputs                | the last expression, can converted to Dict{String} |
| prerequisites            | inputs, outputs       | the last expression                                |
| main (JuliaProgram only) | inputs, outputs       | outputs::Dict{String}                              |
| validate_outputs         | outputs               | the last expression                                |
| wrap_up                  | inputs, outputs       | the last expression                                |

### Example

```julia
prog = JuliaProgram(
    inputs = ["A", "B"],
    outputs = ["OUT"],
    validate_inputs = quote
        @show A
        @show inputs
        A isa Number
    end,
    infer_outputs = quote
        Dict("OUT" => A + B)
    end,
    main = quote
        @show A
        @show B
        OUT = A + B
    end,
    validate_outputs = quote
        @show OUT
        OUT isa Number
    end
)

run(prog; A = 3, B = 5, touch_run_id_file = false)
# (true, Dict{String, Any}("OUT" = 8))
```

!!! warning "`quote` variables in other scopes"
    1. A local variable (include function) should be referenced using `\$` in `expr`ession. (No need to use `\$` for global variables.)

    2. A local `::Symbol` variable (`sym`) should be referenced using `\$(QuoteNode(sym))` in `expr`ession.

    ### Example:
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
            @show \$(QuoteNode(l_sym))
            @show \$l_var + 2
            \$l_func()
            A + B
        end
    end
    
    expr = gen_expr()
    func = Pipelines.quote_function(expr, inputs; mod = @__MODULE__)
    
    in_dict = Dict("A" => 5, "B" => 50)
    func(in_dict)
    ```

!!! compat "Compatibility of Pipelines < v0.8"
    You can still pass `Function` to variables require `Expr`, but you **cannot** use the 'elements as variables' feature. The function should take `inputs::Dict{String}` and/or `outputs::Dict{String}` as variables, and you have to use traditional `inputs["VARNAME"]` in functions.

See also: [`CmdProgram`](@ref), [`JuliaProgram`](@ref), [`quote_function`](@ref)
"""
function quote_expr end

"""

```julia
quote_function(expr::Expr, inputs::Vector{String}; specific_return = nothing, mod::Module = Pipelines)
    # Return `Function` with one argument: `inputs::Dict{String}`.

quote_function(expr::Expr, inputs::Vector{String}, outputs::Vector{String}; specific_return = nothing, mod::Module = Pipelines)
    # Return `Function` with two arguments: `inputs::Dict{String}, outputs::Dict{String}`.

quote_function(f::Function, x; specific_return, mod) = f
quote_function(f::Function, x, y; specific_return, mod) = f
    # Directly return `f::Function` without any process.
```

### Description

When building `Program`, `Expr` are automatically converted to `Function` using `quote_function`. The elements of `inputs` and/or `outputs` in `expr` will be replaced by `inputs["element"]` and/or `outputs["elements"]`, respectively. Also, in the generated function, the arguments (inputs and outputs) are regarded as `Dict{String}`.

- `specific_return`: an `Expr` appended to `expr`.

- `mod::Module`: `Expr`ressions will evaluated to functions in `mod`. Please use `mod = @__MODULE__` to prevent precompilation fail when defining the program within a package.

In `expr::Expr`, elements in `inputs` or `outputs` can be directly used as variables for those arguments. See the table below.

| Argument                 | Elements as variables | Default returned value                             |
| :----------------------- | :-------------------- | :------------------------------------------------- |
| validate_inputs          | inputs                | the last expression                                |
| infer_outputs            | inputs                | the last expression, can be converted to Dict{String} |
| prerequisites            | inputs, outputs       | the last expression                                |
| main (JuliaProgram only) | inputs, outputs       | outputs::Dict{String}                              |
| validate_outputs         | outputs               | the last expression                                |
| wrap_up                  | inputs, outputs       | the last expression                                |

### Usage in Program building

```julia
function JuliaProgram(; kwargs...)
    ...
    # inputs isa Vector{String}
    # outputs isa Vector{String}
    # mod isa Module where evaluating expressions to functions in

    validate_inputs = quote_function(validate_inputs, inputs; mod = mod)
    infer_outputs = quote_function(infer_outputs, inputs; mod = mod)
    prerequisites = quote_function(prerequisites, inputs, outputs; mod = mod)
    validate_outputs = quote_function(validate_outputs, outputs; mod = mod)
    wrap_up = quote_function(wrap_up, inputs, outputs; mod = mod)

    main = quote_function(main, inputs, outputs; specific_return = :(outputs), mod = mod)
    ...
end
```
"""
function quote_function(expr::Expr, inputs::Vector{String}; specific_return = nothing, mod::Module = @__MODULE__)
    expr = dictreplace_all!(expr, inputs, :inputs)
    if isnothing(specific_return)
        @eval mod function(inputs::Dict{String})
            $expr
        end
    else
        specific_return = dictreplace_all!(specific_return, inputs, :inputs)
        @eval mod function(inputs::Dict{String})
            $expr
            $specific_return
        end
    end
end
function quote_function(expr::Expr, inputs::Vector{String}, outputs::Vector{String}; specific_return = nothing, mod::Module = @__MODULE__)
    expr = dictreplace_all!(expr, inputs, :inputs)
    expr = dictreplace_all!(expr, outputs, :outputs)
    if isnothing(specific_return)
        @eval mod function(inputs::Dict{String}, outputs::Dict{String})
            $expr
        end
    else
        specific_return = dictreplace_all!(specific_return, inputs, :inputs)
        specific_return = dictreplace_all!(specific_return, outputs, :outputs)
        @eval mod function(inputs::Dict{String}, outputs::Dict{String})
            $expr
            $specific_return
        end
    end
end
quote_function(f::Function, x::Vector{String}; args...) = f
quote_function(f::Function, x::Vector{String}, y::Vector{String}; args...) = f

dictreplace!(ex, s::Symbol, v::Expr) = ex
dictreplace!(ex::Symbol, s::Symbol, v::Expr) = s == ex ? v : ex
function dictreplace!(ex::AbstractString, s::Symbol, v::Expr)
    if '$' in ex && occursin(string(s), ex)
        ex = Meta.parse(string('"', ex, '"'))  # expose the Expr of $(...)
        ex = dictreplace!(ex, s, v)
        ex = string(ex)[2:end-1]  # has to convert to string, otherwise @eval throws error.
                                  # [2:end-1] removes " at the begin and end
    else
        ex
    end
end
function dictreplace!(ex::Expr, s::Symbol, v::Expr)
    if ex.head == :kw
        for i in 2:length(ex.args)  # func(a; b = b): the first b does not change
            ex.args[i] = dictreplace!(ex.args[i], s, v)
        end
    else
        for i in 1:length(ex.args)
            ex.args[i] = dictreplace!(ex.args[i], s, v)
        end
    end
    ex
end

function dictreplace_all!(expr, kys, dsym::Symbol)
    for k in kys
        expr = dictreplace!(expr, Symbol(k), :($(dsym)[$k]))
    end
    expr
end
