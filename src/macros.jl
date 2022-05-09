"""
    @vars program::Program args...

Generate runable `(inputs::Dict, outputs::Dict)` for `program` using `args` in the form of `key = value`.

## Example
```julia
p = CmdProgram(
	inputs = [
		"input",
		"input2" => Int,
		"optional_arg" => 5,
		"optional_arg2" => 0.5 => Number
	],
	outputs = [
		"out" => "<input>.output"
	],
	cmd = `echo input input2 optional_arg optional_arg2 output`
)

i = "iout"
kk = :xxx
inputs, outputs = @vars p input=kk input2=22 optional_arg=:sym out=i
run(p, inputs, outputs)
```
"""
macro vars(program, args...)
    return quote
        local p = $(esc(program))
        if !(p isa Program)
            error("The first argument of @vars should be a ::Program")
        end
        local invars = Symbol.(p.inputs)
        local outvars = Symbol.(p.outputs)
        local inputs = Dict{String,Any}()
        local outputs = Dict{String,Any}()
        local narg = length($args)
        local i = 1
        while i <= narg
            local arg = ($args)[i]
            if arg.head === :(=)
                local key = arg.args[1]
                local val = arg.args[2]
                if key in invars
                    if val isa QuoteNode # val is Symbol
                        setindex!(inputs, eval(val), string(key))
                    elseif val isa Symbol
                        setindex!(inputs, Core.eval(@__MODULE__, val), string(key))
                    else
                        setindex!(inputs, val, string(key))
                    end
                elseif key in outvars
                    # same as invars
                    if val isa QuoteNode # val is Symbol
                        setindex!(outputs, eval(val), string(key))
                    elseif val isa Symbol
                        setindex!(outputs, Core.eval(@__MODULE__, val), string(key))
                    else
                        setindex!(outputs, val, string(key))
                    end
                end
            end
            i += 1
        end
        inputs, outputs
    end
end
