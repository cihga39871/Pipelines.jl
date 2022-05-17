
## macro
jp = JuliaProgram(
	name = "Echo",
	id_file = "id_file",
	inputs = [
		"input1",
		"input2" => Int,
		"optional_arg" => 5,
		"optional_arg2" => 0.5 => Number
	],
	outputs = [
		"output" => "<input1>.output"
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
prog_run(jp; input1 = kk, input2 = 33, output = i, commonargs...)

## v0.7.6 @extractvars wrapped in p.main(inputs, outputs)
#=jp2 = JuliaProgram(
	name = "Echo2",
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
		@show input
		@show input2
		@show optional_arg
		@show optional_arg2
		@show output
	end
)

main = (inputs, outputs) -> begin
	localvars(inputs, outputs, quote
		@show X
		@show Y
		@show Z
	end)
end


i = Dict("X" => :in, "Y" => 22)
o = Dict("Z" => "pit")

main(i,o)

@show input2


=#
