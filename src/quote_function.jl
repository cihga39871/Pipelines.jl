
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
    1. A local variable in other scopes should be referenced using `\$` in `expr`ession. (No need to use `\$` for global variables.)

    2. A local `::Symbol` variable (`sym`) should be referenced using `\$(QuoteNode(sym))` in `expr`ession.

    ### Example:
    ```julia
    inputs = ["A", "B"]
    g_var = 3
    g_sym = :globalsymbol

    function gen_expr()
        l_var = 5
        l_sym = :abc
        expr = quote
            @show inputs
            @show g_var
            @show g_sym
            @show \$(QuoteNode(l_sym))
            @show \$l_var + 2
            A + B
        end
    end

    expr = gen_expr()
    ```

!!! compat "Compatibility of Pipelines < v0.8"
    You can still pass `Function` to variables require `Expr`, but you **cannot** use the 'elements as variables' feature. The function should take `inputs::Dict{String}` and/or `outputs::Dict{String}` as variables, and you have to use traditional `inputs["VARNAME"]` in functions.

See also: [`CmdProgram`](@ref), [`JuliaProgram`](@ref), [`quote_function`](@ref)
"""
function quote_expr end

"""
    quote_function(expr::Expr, inputs::Vector{String}; specific_return = nothing)

Return `Function` with one argument: `inputs::Dict{String}`

```julia
quote_function(expr::Expr, inputs::Vector{String}, outputs::Vector{String}; specific_return = nothing)
```

Return `Function` with two arguments: `inputs::Dict{String}, outputs::Dict{String}`

```julia
quote_function(f::Function, x; specific_return) = f
quote_function(f::Function, x, y; specific_return) = f
```

Directly return `f::Function` without any process.

### Description

When building `Program`, `Expr` are automatically converted to `Function` using `quote_function`. The elements of `inputs` and/or `outputs` in `expr` will be replaced by `inputs["element"]` and/or `outputs["elements"]`, respectively. Also, in the generated function, the arguments (inputs and outputs) are regarded as `Dict{String}`.

- `specific_return`: an `Expr` appended to `expr`.

In `expr::Expr`, elements in `inputs` or `outputs` can be directly used as variables for those arguments. See the table below.

| Argument                 | Elements as variables | Default returned value                             |
| :----------------------- | :-------------------- | :------------------------------------------------- |
| validate_inputs          | inputs                | the last expression                                |
| infer_outputs            | inputs                | the last expression, can converted to Dict{String} |
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

    validate_inputs = quote_function(validate_inputs, inputs)
    infer_outputs = quote_function(infer_outputs, inputs)
    prerequisites = quote_function(prerequisites, inputs, outputs)
    validate_outputs = quote_function(validate_outputs, outputs)
    wrap_up = quote_function(wrap_up, inputs, outputs)

    main = quote_function(main, inputs, outputs; specific_return = :(outputs))
    ...
end
```

!!! warning "`quote` variables in other scopes"
    1. A local variable in other scopes should be referenced using `\$` in `expr`ession. (No need to use `\$` for global variables.)

    2. A local `::Symbol` variable (`sym`) should be referenced using `\$(QuoteNode(sym))` in `expr`ession.

    ### Example:
    ```julia
    inputs = ["A", "B"]
    g_var = 3
    g_sym = :globalsymbol

    function gen_expr()
        l_var = 5
        l_sym = :abc
        expr = quote
            @show inputs
            @show g_var
            @show g_sym
            @show \$(QuoteNode(l_sym))
            @show \$l_var + 2
            A + B
        end
    end

    expr = gen_expr()
    new_func = quote_function(expr, inputs)

    inputs_of_new_func = Dict("A" => 100, "B" => 50)
    x = new_func(inputs_of_new_func)
    # inputs = Dict("B" => 50, "A" => 100)
    # g_var = 3
    # g_sym = :globalsymbol
    # :abc = :abc
    # 5 + 2 = 7
    # 150

    x
    # 150
    ```
"""
function quote_function(expr::Expr, inputs::Vector{String}; specific_return = nothing, mod::Module = @__MODULE__)
    @show @__MODULE__
    @show parentmodule(@__MODULE__)
    @show mod
    expr = dictreplace_all!(expr, inputs, :inputs)
    if isnothing(specific_return)
        @eval mod function(inputs)
            $expr
        end
    else
        specific_return = dictreplace_all!(specific_return, inputs, :inputs)
        @eval mod function(inputs)
            $expr
            $specific_return
        end
    end
end
function quote_function(expr::Expr, inputs::Vector{String}, outputs::Vector{String}; specific_return = nothing, mod::Module = @__MODULE__)
    @show @__MODULE__
    @show parentmodule(@__MODULE__)
    expr = dictreplace_all!(expr, inputs, :inputs)
    expr = dictreplace_all!(expr, outputs, :outputs)
    if isnothing(specific_return)
        @eval mod function(inputs, outputs)
            $expr
        end
    else
        specific_return = dictreplace_all!(specific_return, inputs, :inputs)
        specific_return = dictreplace_all!(specific_return, outputs, :outputs)
        @eval mod function(inputs, outputs)
            $expr
            $specific_return
        end
    end
end
quote_function(f::Function, x::Vector{String}; specific_return = nothing) = f
quote_function(f::Function, x::Vector{String}, y::Vector{String}; specific_return = nothing) = f

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
    for i in 1:length(ex.args)
        ex.args[i] = dictreplace!(ex.args[i], s, v)
    end
    ex
end

function dictreplace_all!(expr, kys, dsym::Symbol)
    for k in kys
        expr = dictreplace!(expr, Symbol(k), :($(dsym)[$k]))
    end
    expr
end
