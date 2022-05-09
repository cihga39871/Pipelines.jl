"""
    @vars program::Program key_value_args...

Return runable `(inputs::Dict, outputs::Dict)` for `program` using `key_value_args` in the form of `key = value`.

See also [`@run`](@ref)

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
run(p, inputs, outputs; touch_run_id_file = false)
# (true, Dict{String, Any}("output" => "iout"))
```
"""
macro vars(program, args...)
    return quote
        local p = $(esc(program))
        if !(p isa Program)
            error("The first argument of @vars should be a ::Program")
        end
        local inputs = Dict{String,Any}()
        local outputs = Dict{String,Any}()
        local narg = length($args)
        local i = 1
        while i <= narg
            local arg = ($args)[i]
            if arg.head === :(=)
                local key = string(arg.args[1])
                local val = arg.args[2]
                if key in p.inputs
                    setindex!(inputs, Core.eval(@__MODULE__, val), string(key))
                elseif key in p.outputs
                    setindex!(outputs, Core.eval(@__MODULE__, val), string(key))
                else  # may be keyword parameters of other functions, such as run, Job.
                    error("Program $(p.name) does not have a key named $key")
                end
            else
                error("Syntax error: args only support `key = value` form: $arg")
            end
            i += 1
        end
        inputs, outputs
    end
end

"""
    @run program::Program key_value_args... run_args...

Run `program` without creating `inputs::Dict` and `outputs::Dict`.

- `key_value_args`: the inputs and outputs are provided in the form of `key = value`, rather than `Dict`.

- `run_args`: the keyword arguments pass to `run(p::Program, inputs, outputs, run_args...)`.

See also [`@vars`](@ref) [`run`](@ref)

### Example
```julia
jp = JuliaProgram(
	name = "Echo",
	id_file = "id_file",
	inputs = [
		"input",
		"input2" => Int,
		"optional_arg" => 5,
		"optional_arg2" => 0.5 => Number
	],
	outputs = [
		"output" => "<input>.output"
	],
	main = (x,y) -> begin
		@show x
		@show y
		y
	end
)

i = "iout"
kk = :xxx
b = false
commonargs = (touch_run_id_file = b, verbose = :min)
run_res = @run jp input=kk input2=22 optional_arg=:min output=i touch_run_id_file=b
run_res = @run jp input=kk input2=22 optional_arg=:min output=i touch_run_id_file=b verbose=:min
run_res = @run jp input=kk input2=22 optional_arg=:sym output=i commonargs...
# (true, Dict{String, Any}("output" => "iout"))
```
"""
macro run(program, args...)
    return quote
        local p = $(esc(program))
        if !(p isa Program)
            error("The first argument of @vars should be a ::Program")
        end
        local inputs = Dict{String,Any}()
        local outputs = Dict{String,Any}()
        local args = [($args)...]
        local narg = length(args)
        local kw_args = Vector{Any}()  # keyword parameters of other functions, such as run, Job.
        local i = 1
        while i <= narg
            local arg = args[i]
            if arg.head === :(=)
                local key = string(arg.args[1])
                local val = arg.args[2]

                if key in p.inputs
                    setindex!(inputs, Core.eval(@__MODULE__, val), key)
                elseif key in p.outputs
                    setindex!(outputs, Core.eval(@__MODULE__, val), key)
                else  # may be keyword parameters of other functions, such as run, Job.
                    arg.head = :kw  # keyword head is :kw, rather than :(=)
                    if !(val isa QuoteNode) # QuoteNode: such as :(:sym)
                        arg.args[2] = Core.eval(@__MODULE__, val)
                    end
                    push!(kw_args, arg)
                end
            elseif arg.head == :(...)  # common_args...
                local extra_args = Core.eval(@__MODULE__, arg.args[1])
                local extra_args_keys = keys(extra_args)
                local extra_args_vals = values(extra_args)
                local nextra = length(extra_args_keys)
                for m in 1:nextra
                    local k = extra_args_keys[m]
                    local v = extra_args_vals[m]
                    if v isa Symbol
                        v = QuoteNode(v)  # no interpolation of Symbol
                    end
                    push!(args, Expr(:(=), k, v))
                end
                narg = length(args)
            else
                error("SyntaxError: args only support `key = value` form: $arg")
            end
            i += 1
        end

        local ex_kw = Expr(:parameters, kw_args...)
        local ex = :(run($p, $inputs, $outputs))
        insert!(ex.args, 2, ex_kw)  # insert keyword arguments to ex
        eval(ex)
    end
end
