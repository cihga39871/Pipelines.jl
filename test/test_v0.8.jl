
@testset "v0.8" begin
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
	@test run(jp; input1 = kk, input2 = 33, output = i, commonargs...) == (true, Dict{String, Any}("output" => "iout"))

	@test run(jp, Dict(), "output" => "old"; input1 = kk, input2 = 33, output = i, commonargs...) == (true, Dict{String, Any}("output" => "iout"))

	## v0.8.0 @extractvars wrapped in p.main(inputs, outputs)

	# example of CmdProgram
	program_bowtie2 = CmdProgram(
	    name = "Bowtie2 Mapping",
	    id_file = ".bowtie2",

	    inputs = [
	        "FASTQ" => String,
	        "REF" => "human_genome_hg38.fa" => String
	    ],

		outputs = ["BAM" => String],
		infer_outputs = quote
	        Dict("BAM" => to_str(FASTQ) * ".bam")
	    end,

		validate_inputs = quote
	        check_dependency_file(FASTQ) && check_dependency_file(REF)
	    end,

		prerequisites = quote
	        mkpath(dirname(to_str(BAM)))
	    end,

		validate_outputs = quote
	        check_dependency_file(BAM)
	    end,

	    cmd = pipeline(`bowtie2 -x REF -q FASTQ`, `samtools sort -O bam -o BAM`),

	    wrap_up = quote
			run(`samtools index $BAM`)
		end
	)

	p = CmdProgram(
	  id_file = "id_file",
	  inputs = ["input",
				"input2" => Int,
				"optional_arg" => 5,
				"optional_arg2" => 0.5 => Number],
	  outputs = "output" => "<input>.output",
	  validate_inputs = quote
		  optional_arg2 isa Float64 && inputs isa Dict
	  end,
	  cmd = `echo input input2 optional_arg optional_arg2 output`)

	# running the program: keyword arguments include keys of inputs and outputs
	success, outputs = run(p; input = `in1`, input2 = 2, output = "out", touch_run_id_file = false)
	@test success


	# example of JuliaProgram
	p = JuliaProgram(
	  	id_file = "id_file",
	  	inputs = ["a",
	  	          "b" => Int],
	  	outputs = "c" => "<a>.<b>",
	  	main = quote
	  		println("inputs are ", a, " and ", b)
	  		println("You can also use info in outputs: ", c)
	          println("The returned value will be assigned to a new outputs")
	  		c = b^2
	  	end)

	# running the program: keyword arguments include keys of inputs and outputs
	success, new_out = run(p; a = `in1`, b = 2, c = "out", touch_run_id_file = false)
	@test success

	@test new_out != infer_outputs(p; a = `in1`, b = 2, c = "out")  # outputs will change to the returned value of main function, if the returned value is a Dict and pass `p.validate_outputs`


	# example of quote_expr
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

	@test run(prog; A = 3, B = 5, touch_run_id_file = false)[1]
	# (true, Dict{String, Any}("OUT" = 8))

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
			@show $(QuoteNode(l_sym))
			@show $l_var + 2
			A + B
		end
	end

	@test_nowarn expr = gen_expr()

end
