
"""
    quote_function(expr::Expr, inputs::Vector{String})
    quote_function(expr::Expr, inputs::Vector{String}, outputs::Vector{String})

Generate function using `expr`. The elements of `inputs` and `outputs` in `expr` will be replaced by `inputs["element"]`. Also, in `expr`, inputs and outputs are regarded as a `Dict{String}`.

Variables from other scopes used in `expr` should be referenced using `\$`, such as:

    x = 5
    expr = quote \$x + 2 end

### Example

```julia

    inputs = ["A", "B"]
    var_in_other_scope = "use_dollar_sign"
	symbol = :symbol_need_QuoteNode
    expr = quote
        @show A + B
        @show inputs
        @show $var_in_other_scope
		@show $(QuoteNode(symbol))
		return "returned_value"
    end

    new_func = quote_function(expr, inputs)

    new_func(Dict("A" => 3, "B" => 33))
    # inputs["A"] + inputs["B"]
    # inputs = Dict("A" => 3, "B" => 33)
    # "use_dollar_sign" = "use_dollar_sign"
	# :symbol_need_QuoteNode = :symbol_need_QuoteNode
```
"""
function quote_function(expr::Expr, inputs::Vector{String})
    expr_new = dictreplace_all(expr, inputs, :inputs)
	@eval function(inputs)
		$expr_new
	end
end


dictreplace!(ex, s, v) = ex
dictreplace!(ex::Symbol, s, v) = s == ex ? v : ex
function dictreplace!(ex::Expr, s, v)
	for i in 1:length(ex.args)
		ex.args[i] = dictreplace!(ex.args[i], s, v)
	end
	ex
end
dictreplace(ex, s, v) = dictreplace!(copy(ex), s, v)

function dictreplace_all(expr, kys, dsym)
    for k in kys
        expr = dictreplace(expr, Symbol(k), :($(dsym)[$k]))
    end
	expr
end
