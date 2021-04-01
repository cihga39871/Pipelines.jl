var documenterSearchIndex = {"docs":
[{"location":"command_dependency/#Command-Dependency","page":"Command Dependency","title":"Command Dependency","text":"","category":"section"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"Command Dependency (CmdDependency) is the external dependent program that can be used in Command Program (CmdProgram).","category":"page"},{"location":"command_dependency/#Create","page":"Command Dependency","title":"Create","text":"","category":"section"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"CmdDependency(;\n    exec::Base.AbstractCmd=``,\n    test_args::Base.AbstractCmd=``,\n    validate_success::Bool=false,\n    validate_stdout::Function=do_nothing,\n    validate_stderr::Function=do_nothing,\n    exit_when_fail::Bool=true\n)","category":"page"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"exec::AbstractCmd: the command to call the dependency.\ntest_args::AbstractCmd: for testing purposes, the command to be appended to exec.\nvalidate_success::Bool: when checking the dependency, whether to validate the exit code == 0.\nvalidate_stdout::Function: a function takes standard out as String and return the validation result as ::Bool.\nvalidate_stderr::Function: a function takes standard error as String and return the validation result as ::Bool.\nexit_when_fail::Bool: if validation fails, whether to throw error and exit.","category":"page"},{"location":"command_dependency/#Example","page":"Command Dependency","title":"Example","text":"","category":"section"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"julia = CmdDependency(\n    exec = `julia`,\n    test_args = `--version`,\n    validate_success = true,\n    validate_stdout = x -> occursin(r\"^julia version\", x),\n    validate_stderr = do_nothing,\n    exit_when_fail = true\n)","category":"page"},{"location":"command_dependency/#Check","page":"Command Dependency","title":"Check","text":"","category":"section"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"CmdDependency can be checked manually by using","category":"page"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"check_dependency(p::CmdDependency)","category":"page"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"Usually, it is not necessary when you specify CmdDependency in CmdProgram because dependency check will be done when running CmdProgram by default.","category":"page"},{"location":"command_dependency/#Use-with-CmdProgram","page":"Command Dependency","title":"Use with CmdProgram","text":"","category":"section"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"When building CmdProgram, we can specify our command dependencies by using","category":"page"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"dep1 = CmdDependency(...)\ndep2 = CmdDependency(...)\n\nprog = CmdProgram(\n\tcmd_dependencies = [dep1, dep2],\n    ...\n)","category":"page"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"To call the dependency in command, we can use exec(dep) or dep.exec:","category":"page"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"prog = CmdProgram(\n\tcmd_dependencies = [dep1, dep2],\n\tcmd = `$(exec(dep1)) --args` & `$(dep2.exec) --args`\n)","category":"page"},{"location":"command_dependency/#An-Example","page":"Command Dependency","title":"An Example","text":"","category":"section"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"We use the same example in Command Program page.","category":"page"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"Adding samtools and bowtie2 as the dependencies of the bowtie2 mapping program:","category":"page"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"SAMTOOLS = CmdDependency(\n\texec = `samtools`,\n\ttest_args = `--version`,\n\tvalidate_success = true,\n\tvalidate_stdout = x -> occursin(r\"^samtools \\d+\", x)\n)\n\nBOWTIE2 = CmdDependency(\n\texec = `bowtie2`,\n\ttest_args = `--version`,\n\tvalidate_success = true,\n\tvalidate_stdout = x -> occursin(r\"bowtie2-align-s version \\d+\", x)\n)\n\nprogram_bowtie2 = CmdProgram(\n\tname = \"Bowtie2 Mapping\",\n\tid_file = \".bowtie2\",\n\n\tcmd_dependencies = [SAMTOOLS, BOWTIE2],\n\n\tinputs = [\"FASTQ\", \"REF\"],\n\tvalidate_inputs = inputs -> begin\n\t\tcheck_dependency_file(inputs[\"FASTQ\"]) &&\n\t\tcheck_dependency_file(inputs[\"REF\"])\n\tend,\n\n\tprerequisites = (inputs, outputs) -> begin\n\t\tmkpath(dirname(to_str(outputs[\"BAM\"])))\n\tend,\n\n\toutputs = [\"BAM\"],\n\tinfer_outputs = inputs -> begin\n\t\tDict(\"BAM\" => str(inputs[\"FASTQ\"]) * \".bam\")\n\tend,\n\tvalidate_outputs = outputs -> begin\n\t\tcheck_dependency_file(outputs[\"BAM\"])\n\tend,\n\n\tcmd = pipeline(`$(BOWTIE2.exec) -x REF -q FASTQ`, `$(SAMTOOLS.exec) sort -O bam -o BAM`),\n\n\twrap_up = (inputs, outputs) -> run(`$(SAMTOOLS.exec) index $(outputs[\"BAM\"])`)\n)","category":"page"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"To run the program, we can decide whether to check dependencies with run(..., check_dependencies=true):","category":"page"},{"location":"command_dependency/","page":"Command Dependency","title":"Command Dependency","text":"inputs = Dict(\"FASTQ\" => \"a.fastq\", \"REF\" => \"ref.fasta\")\nsuccess, outputs = run(program_bowtie2, inputs; check_dependencies = true)","category":"page"},{"location":"command_program/#Command-Program","page":"Command Program","title":"Command Program","text":"","category":"section"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"Please read the Quick Start section in Home page first.","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"Pipelines are built with multiple CmdPrograms. A CmdProgram contains a command template, lists of dependencies/inputs/outputs, and predefined methods to prepare run-time environment, validate dependencies/inputs/outputs, and so on.","category":"page"},{"location":"command_program/#An-Example","page":"Command Program","title":"An Example","text":"","category":"section"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"We will go through an example to illustrate how to write a CmdProgram.","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"The example is a robust Bowtie2 mapping program. It shows every functionality of CmdProgram, but in reality, we do not need that much validation.","category":"page"},{"location":"command_program/#Define-the-Name","page":"Command Program","title":"Define the Name","text":"","category":"section"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"program_bowtie2 = CmdProgram(\n\tname = \"Bowtie2 Mapping\",# The name of CmdProgram\n\tid_file = \".bowtie2\"     # When job completed, a file \".bowtie2.xxxxxx\" will\n\t                         # be created to indicate the job is finished to\n\t\t\t\t\t\t\t # avoid re-run.\n)","category":"page"},{"location":"command_program/#Main-Command-(Required)","page":"Command Program","title":"Main Command (Required)","text":"","category":"section"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"In Bash script, the main code is","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"REF=/path/to/reference_genome\nFASTQ=/path/to/input_fastq_file\nBAM=/path/to/output_bam_file\n\nbowtie2 -x $REF -q $FASTQ | samtools sort -O bam -o $BAM","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"Code in CmdProgram:","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"program_bowtie2 = CmdProgram(\n\t...,\n\tinputs = [\"FASTQ\", \"REF\"],\n\toutputs = [\"BAM\"],\n\tcmd = pipeline(`bowtie2 -x REF -q FASTQ`, `samtools sort -O bam -o BAM`)\n)","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"Now, the code can be run by invoking run(program_bowtie2, inputs, outputs), but to illustrate the full functionality, we will add more things to make it robust and easy to use.","category":"page"},{"location":"command_program/#Command-Dependency-(Robustness)","page":"Command Program","title":"Command Dependency (Robustness↑)","text":"","category":"section"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"We use samtools and bowtie2 as command dependencies. They can be wrapped in CmdDependency, which is illustrated in another page. Please allow me to skip it for now.","category":"page"},{"location":"command_program/#Infer-Outputs-(Convenience)","page":"Command Program","title":"Infer Outputs (Convenience↑)","text":"","category":"section"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"We can set a default method to generate outputs::Dict{String} from inputs, which allows us run the program without specifying outputs (run(program_bowtie2, inputs).)","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"program_bowtie2 = CmdProgram(\n\t...,\n\toutputs = [\"BAM\"],\n\tinfer_outputs = inputs -> begin\n\t\tDict(\"BAM\" => to_str(inputs[\"FASTQ\"]) * \".bam\")\n\tend,\n\t...\n)","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"note: to_str(x) and to_cmd(x)\nto_str converts most types to String, and to_cmd to Cmd. They are tailored for parsing inputs[\"x\"] and outputs[\"x\"].User-defined inputs, outputs::Dict{String} only confine the key type (String), does not confine the value type because of flexibility. When writing functions using inputs/outputs, we should consider this. It can be a Number, a Cmd, a String, and even a Vector. Pipeline.jl provides to_str and to_cmd to elegantly convert those types to String or Cmd as you wish.Other conversions are also available, such as replaceext (replace extension) and removeext (remove extension). More details are in API/Utils page.","category":"page"},{"location":"command_program/#Validate-Inputs-(Robustness)","page":"Command Program","title":"Validate Inputs (Robustness↑)","text":"","category":"section"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"To make the code robust, we can check whether the inputs exists by using validate_inputs. It is a function takes inputs::Dict{String} as the argument.","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"program_bowtie2 = CmdProgram(\n\t...,\n\tinputs = [\"FASTQ\", \"REF\"],\n\tvalidate_inputs = inputs -> begin\n\t\tcheck_dependency_file(inputs[\"FASTQ\"]) &&\n\t\tcheck_dependency_file(inputs[\"REF\"])\n\tend,\n\t...\n)","category":"page"},{"location":"command_program/#Prerequisites-(Robustness)","page":"Command Program","title":"Prerequisites (Robustness↑)","text":"","category":"section"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"We also need to prepare something (prerequisites) before running the main command. For example, create the output directory if not exist. (prerequisites) is a function takes inputs::Dict{String}, outputs::Dict{String} as the arguments.","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"program_bowtie2 = CmdProgram(\n\t...,\n\tprerequisites = (inputs, outputs) -> begin\n\t\tmkpath(dirname(to_str(outputs[\"BAM\"])))\n\tend\n)","category":"page"},{"location":"command_program/#Validate-Outputs-(Robustness)","page":"Command Program","title":"Validate Outputs (Robustness↑)","text":"","category":"section"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"After running the main command, we can validate outputs by using validate_outputs. It is a function takes outputs::Dict{String} as the argument.","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"program_bowtie2 = CmdProgram(\n\t...,\n\toutputs = [\"BAM\"],\n\tvalidate_outputs = outputs -> begin\n\t\tcheck_dependency_file(outputs[\"BAM\"])\n\tend,\n\t...\n)","category":"page"},{"location":"command_program/#Wrap-Up-(Convenience)","page":"Command Program","title":"Wrap Up (Convenience↑)","text":"","category":"section"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"After validating outputs, we may also do something to wrap up, such as removing temporary files. Here, we build an index for output BAM file. wrap_up function takes inputs::Dict{String}, outputs::Dict{String} as the arguments.","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"program_bowtie2 = CmdProgram(\n\t...,\n\twrap_up = (inputs, outputs) -> run(`samtools index $(outputs[\"BAM\"])`)\n)","category":"page"},{"location":"command_program/#The-Final-Code","page":"Command Program","title":"The Final Code","text":"","category":"section"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"All in all, the final program is like this:","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"program_bowtie2 = CmdProgram(\n\tname = \"Bowtie2 Mapping\",\n\tid_file = \".bowtie2\",\n\n\tinputs = [\"FASTQ\", \"REF\"],\n\tvalidate_inputs = inputs -> begin\n\t\tcheck_dependency_file(inputs[\"FASTQ\"]) &&\n\t\tcheck_dependency_file(inputs[\"REF\"])\n\tend,\n\n\tprerequisites = (inputs, outputs) -> begin\n\t\tmkpath(dirname(to_str(outputs[\"BAM\"])))\n\tend,\n\n\toutputs = [\"BAM\"],\n\tinfer_outputs = inputs -> begin\n\t\tDict(\"BAM\" => str(inputs[\"FASTQ\"]) * \".bam\")\n\tend,\n\tvalidate_outputs = outputs -> begin\n\t\tcheck_dependency_file(outputs[\"BAM\"])\n\tend,\n\n\tcmd = pipeline(`bowtie2 -x REF -q FASTQ`, `samtools sort -O bam -o BAM`),\n\n\twrap_up = (inputs, outputs) -> run(`samtools index $(outputs[\"BAM\"])`)\n)","category":"page"},{"location":"command_program/#Structure","page":"Command Program","title":"Structure","text":"","category":"section"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"CmdProgram can be build with the method:","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"CmdProgram(;\n\tname::String               = \"Unnamed Command Program\",\n\tid_file::String            = \"\",\n\tinfo_before::String        = \"auto\",\n\tinfo_after::String         = \"auto\",\n\tcmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),\n\tinputs::Vector{String}     = Vector{String}(),\n\tvalidate_inputs::Function  = do_nothing,  # positional arguments: inputs::Dict{String}\n\tprerequisites::Function    = do_nothing,  # positional arguments: inputs, outputs::Dict{String}\n\tcmd::Base.AbstractCmd      = ``,\n\tinfer_outputs::Function    = do_nothing,  # positional arguments: inputs::Dict{String}\n\toutputs::Vector{String}    = Vector{String}(),\n\tvalidate_outputs::Function = do_nothing,  # positional arguments: outputs::Dict{String},\n\twrap_up::Function          = do_nothing   # positional arguments: inputs, outputs::Dict{String}\n) -> CmdProgram","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"In this way, all preparation and post-evaluation can be wrapped in a single CmdProgram. It is easy to maintain and use.","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"To run a CmdProgram, use one of the following methods:","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"run(\n\tp::CmdProgram;\n\tinputs::Dict{String}=Dict{String, Cmd}(),\n\toutputs::Dict{String}=Dict{String, Cmd}(),\n\tskip_when_done::Bool=true,\n\tcheck_dependencies::Bool=true,\n\tstdout=nothing,\n\tstderr=nothing,\n\tappend::Bool=false,\n\tverbose::Bool=true,\n\ttouch_run_id_file::Bool=true,\n\tdry_run::Bool=false\n) -> (success::Bool, outputs::Dict{String})\n\nrun(\n\tp::CmdProgram,\n\tinputs::Dict{String},\n\toutputs::Dict{String};\n\tkwargs...\n)\n\nrun(\n\tp::CmdProgram,\n\tinputs::Dict{String},\n\tkwargs...\n)  # only usable when `p.infer_outputs` is defined.","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"The explanation of arguments is in the next section.","category":"page"},{"location":"command_program/#Workflow","page":"Command Program","title":"Workflow","text":"","category":"section"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"Check keywords consistency: Inputs/outputs keywords should be consistent in both p::CmdProgram and run(p; inputs, outputs).\nFor example, if inputs and outputs in p::CmdProgram is defined this wayp = CmdProgram(..., inputs = [\"I\", \"J\"], outputs = [\"K\"])The inputs and outputs in run(...) has to berun(p;\n    inputs = Dict(\n        \"I\" => something,\n        \"J\" => something\n    ),\n    outputs = Dict(\"K\" => something)\n)\nPrint info about starting program.\nThe content can set by p.info_before::String.Disable: run(..., verbose=false)\nCheck whether the program ran successfully before. If so, return (true, outputs::Dict{String}) without running it.\nHow does the program know it ran before?a. run(..., skip_when_down=true) enables the functionality.b. Run id file is generated. The prefix of run id file is set by p.id_file. After given inputs and outputs, a unique ID will be appended to the prefix. You can use run(..., touch_run_id_file=false) to skip creating the run id file.c. p.validate_outputs(outputs) run successfully without returning false.\nCheck command dependencies (CmdDependency).\nDisable: run(..., check_dependencies=false)Read Command Dependency portion for details.\nRemove the run id file if exists.\nValidate inputs by invoking p.validate_inputs(inputs).\nPreparing the main command.\nIf you specify run(...; stdout=something, stderr=something, append::Bool), the command (cmd) will be wrapped with pipeline(cmd; stdout=something, stderr=something, append::Bool). If cmd has its own file redirection, the outer wrapper may not work as you expect.\nMeet prerequisites by invoking p.prerequisites(inputs, outputs).\nIt is the last function before running the main command.  For example, you can create a directory if  the main command cannot create itself.\nRun the main command.\nValidate outputs by invoking p.validate_outputs(outputs).\nRun the wrap up function by invoking p.wrap_up(inputs, outputs)\nIt is the last function to do after-command jobs. For example, you can delete intermediate files if necessary.\nCreate run id file if run(..., touch_run_id_file=true). Read Step 3 for details.\nPrint info about finishing program.\nThe content can set by p.info_after::String.Disable: run(..., verbose=false)\nReturn (success::Bool, outputs{String})","category":"page"},{"location":"command_program/","page":"Command Program","title":"Command Program","text":"note: Dry Run\nrun(..., dry_run=true) will return (mature_command::AbstractCmd, run_id_file::String) instead.","category":"page"},{"location":"#Pipelines.jl","page":"Home","title":"Pipelines.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Pipelines","category":"page"},{"location":"","page":"Home","title":"Home","text":"A lightweight Julia package for computational pipelines.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Building reusable pipelines and workflows is easier than you have ever thought.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Easy to build both simple and complex tasks.\nSupports external command lines and pure Julia functions.\nSupports resuming interrupted tasks, skipping finished tasks.\nSupports dependency check.\nSupports inputs, outputs validation, and so on.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pipelines.jl can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add Pipelines\n# If it fails, use\npkg> add https://github.com/cihga39871/Pipelines.jl.git","category":"page"},{"location":"","page":"Home","title":"Home","text":"To use the package, type","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pipelines","category":"page"},{"location":"#Quick-Start","page":"Home","title":"Quick Start","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pipelines are built with multiple CmdPrograms. A CmdProgram contains a command template and name lists of inputs/outputs. The names of inputs/outputs will be replaced by real values when executing the program.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Let's set up a simple program to print values using echo:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pipelines\n\necho = CmdProgram(\n    inputs = [\"INPUT1\", \"INPUT2\"],\n    cmd = `echo INPUT1 INPUT2`   \n)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Running the program is just like running other Cmd,  but here we need to specify inputs by using Dict{String}.","category":"page"},{"location":"","page":"Home","title":"Home","text":"inputs = Dict(\n    \"INPUT1\" => \"Hello,\",\n    \"INPUT2\" => `Pipeline.jl`\n)\nrun(echo, inputs)","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Program will not run twice by default!\nIf you run a program with the same inputs again, the program will just return the same result, display a warning message without running the command twice.run(echo, inputs)This is because the program will generate a file (run id file) in the current directory indicating the program has been run. Several methods can be used to re-run a program:# Method 1: stop checking finished program\nrun(echo, inputs; skip_when_done = false)\n\n# Method 2: delete the run_id_file before running again:\ncmd, run_id_file = run(echo, inputs; dry_run = true) # Dry-run returns the command and run id file without running it.\nrm(run_id_file)  # remove the run_id_file\n\n# Method 3: Do not generate run_id_file when first running.\nrun(echo, inputs; touch_run_id_file=false)","category":"page"},{"location":"#Program-with-outputs","page":"Home","title":"Program with outputs","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Unlike the first example, many programs write files as outputs. Pipelines.jl has an elegant way to handle it.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The following program prints values simultaneously, sort them, and save to a file.","category":"page"},{"location":"","page":"Home","title":"Home","text":"prog = CmdProgram(\n    inputs = [\"INPUT1\", \"INPUT2\", \"INPUT3\"],\n    outputs = [\"OUTPUT_FILE\"],\n    cmd = pipeline(`echo INPUT1 INPUT2` & `echo INPUT3`, `sort`, \"OUTPUT_FILE\")\n)\n\ninputs = Dict(\n    \"INPUT1\" => \"Hello,\",\n    \"INPUT2\" => `Pipeline.jl`,\n    \"INPUT3\" => 39871\n)\noutputs = Dict(\"OUTPUT_FILE\" => \"out.txt\")\n\nrun(prog, inputs, outputs) # will return (success::Bool, outputs)","category":"page"},{"location":"","page":"Home","title":"Home","text":"It is inconvenient to specify outputs every time, so we provide an argument (infer_outputs::Function) in CmdProgram to generate default outputs from inputs.","category":"page"},{"location":"","page":"Home","title":"Home","text":"prog = CmdProgram(\n    inputs = [\"INPUT1\", \"INPUT2\", \"INPUT3\"],\n    outputs = [\"OUTPUT_FILE\"],\n    cmd = pipeline(`echo INPUT1 INPUT2` & `echo INPUT3`, `sort`, \"OUTPUT_FILE\"),\n    infer_outputs = inputs -> Dict(\n    \t\"OUTPUT_FILE\" => inputs[\"INPUT1\"] * \".txt\"\n    )\n)\nsuccess, outputs = run(prog, inputs)","category":"page"},{"location":"","page":"Home","title":"Home","text":"We can also generate default outputs without running the program:","category":"page"},{"location":"","page":"Home","title":"Home","text":"outputs = infer_outputs(prog, inputs)","category":"page"},{"location":"#Future-development","page":"Home","title":"Future development","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Support running competitive tasks with locks.\nSupport running pure Julia Program.","category":"page"},{"location":"API/#API","page":"API","title":"API","text":"","category":"section"},{"location":"API/#Command-Program","page":"API","title":"Command Program","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"CmdProgram()\nrun(::CmdProgram)\ninfer_outputs(::CmdProgram, ::Dict{String})","category":"page"},{"location":"API/#Pipelines.CmdProgram-Tuple{}","page":"API","title":"Pipelines.CmdProgram","text":"Struct\n\nmutable struct CmdProgram\n\tname::String\n\tid_file::String\n\tinfo_before::String\n\tinfo_after::String\n\tcmd_dependencies::Vector{CmdDependency}\n\tinputs::Vector{String}\n\tvalidate_inputs::Function\n\tprerequisites::Function\n\tcmd::Base.AbstractCmd\n\tinfer_outputs::Function\n\toutputs::Vector{String}\n\tvalidate_outputs::Function\n\twrap_up::Function\nend\n\nMethods\n\nCmdProgram(;\n\tname::String               = \"Unnamed Command Program\",\n\tid_file::String            = \"\",\n\tinfo_before::String        = \"auto\",\n\tinfo_after::String         = \"auto\",\n\tcmd_dependencies::Vector{CmdDependency} = Vector{CmdDependency}(),\n\tinputs::Vector{String}     = Vector{String}(),\n\tvalidate_inputs::Function  = do_nothing,  # positional arguments: inputs::Dict{String, ValidInputTypes}\n\tprerequisites::Function    = do_nothing,  # positional arguments: inputs, outputs::Dict{String, ValidInputTypes}\n\tcmd::Base.AbstractCmd      = ``,\n\tinfer_outputs::Function    = do_nothing,  # positional arguments: inputs::Dict{String, ValidInputTypes}\n\toutputs::Vector{String}    = Vector{String}(),\n\tvalidate_outputs::Function = do_nothing  # positional arguments: outputs::Dict{String, ValidInputTypes},\n\twrap_up::Function          = do_nothing  # positional arguments: inputs, outputs::Dict{String, ValidInputTypes}\n) -> CmdProgram\n\nCommand program template. To run a CmdProgram, use run(::CmdProgram; kwargs...).\n\nArguments\n\nname::String: Program name.\nid_file::String: The prefix of run ID file. To prevent from running the program with the same inputs and outputs twice, it will generate a unique run ID file after a successful run.\ninfo_before::String: Print it when starting the program.\ninfo_after::String: Print it when finishing the program.\ncmd_dependencies::Vector{CmdDependency}: Any command dependencies used in the program.\ninputs and outputs: keywords (Vector{String}) in cmd that can be replaced when envoking run(::CmdProgram, inputs::Dict{String, ValidInputTypes}, outputs::Dict{String, ValidInputTypes}). See details below.\nCmdProgram stores a command template. In the template, replaceable portions are occupied by keywords, and all keywords are set in inputs::Vector{String} and outputs::Vector{String}.\nTo run the program with replaced keywords, you need to use run(::CmdProgram; inputs::Dict{String, ValidInputTypes}, outputs::Dict{String, ValidInputTypes}). The data type is different.\nvalidate_inputs::Function: A function to validate inputs. It takes one argument Dict{String, ValidInputTypes} whose keys are the same as inputs. If validation fail, throw error or return false.\nprerequisites: A function to run just before the main command. It prepares necessary things, such as creating directories. It takes two arguments Dict{String, ValidInputTypes} whose keys are the same as inputs and outputs, respectively.\ncmd::AbstractCmd: The main command template. In the template, keywords in inputs::Vector{String} and outputs::Vector{String} will be replaced when envoking run(::CmdProgram, inputs::Dict{String, ValidInputTypes}, outputs::Dict{String, ValidInputTypes}).\ninfer_outputs::Function: A function to infer outputs from inputs. It takes one argument Dict{String, ValidInputTypes} whose keys are the same as inputs.\nvalidate_outputs::Function: A function to validate outputs. It takes one argument Dict{String, ValidInputTypes} whose keys are the same as outputs. If validation fail, throw error or return false.\nwrap_up::Function: the last function to run. It takes two arguments Dict{String, ValidInputTypes} whose keys are the same as inputs and outputs, respectively.\n\nExample\n\np = CmdProgram(\n\tcmd_dependencies = [julia],\n\tid_file = \"id_file\",\n\tinputs = [\"input\", \"input2\"],\n\toutputs = [\"output\"],\n\tcmd = `echo input input2 output`\n)\n\ninputs = Dict(\n\t\"input\" => `in1`,\n\t\"input2\" => 2\n)\n\noutputs = Dict(\n\t\"output\" => \"out\"\n)\n\nrun(p, inputs, outputs;\n\ttouch_run_id_file = false\n)\n\n\n\n\n\n","category":"method"},{"location":"API/#Base.run-Tuple{CmdProgram}","page":"API","title":"Base.run","text":"run(\n\tp::CmdProgram;\n\tinputs::Dict{String}=Dict{String, Cmd}(),\n\toutputs::Dict{String}=Dict{String, Cmd}(),\n\tskip_when_done::Bool=true,\n\tcheck_dependencies::Bool=true,\n\tstdout=nothing,\n\tstderr=nothing,\n\tappend::Bool=false,\n\tverbose::Bool=true,\n\ttouch_run_id_file::Bool=true,\n\tdry_run::Bool=false\n) -> (success::Bool, outputs::Dict{String})\n\nrun(\n\tp::CmdProgram,\n\tinputs::Dict{String},\n\toutputs::Dict{String};\n\tkwargs...\n)\n\nrun(\n\tp::CmdProgram,\n\tinputs::Dict{String},\n\tkwargs...\n)  # only usable when `p.infer_outputs` is defined.\n\nRun Command Program (CmdProgram) using specified inputs and outputs.\n\nReturn (success::Bool, outputs::Dict{String})\n\np::CmdProgram: the command program template.\ninputs::Dict{String} and outputs::Dict{String}: p::CmdProgram stores a command template. In the template, replaceable portions are occupied by keywords, and all keywords can be found at p.inputs and p.outputs string vectors. Here, inputs and outputs are Dict(keyword::String => replacement::Union{Cmd, AbstractString, Number, Array{T,1} where T}). The replacements do not have a length limit, unless a keyword refers to a filename (length == 1).\nskip_when_done::Bool = true: Skip running the program and return true if it has been done before (the run_id_file exists and p.validate_outputs(outputs) passes.)\ncheck_dependencies::Bool = true: check dependencies for p (p.cmd_dependencies).\nstdout, stderr and append: Redirect the program output to files, the same behavior as pipeline(cmd; stdout=stdout, stderr=stderr, append=append). Caution: use after checking whether the original command has redirection.\nverbose::Bool = true: If true, print info and error messages. If false, print error messages only.\ntouch_run_id_file::Bool = true: If true, touch a unique run ID file, which indicate the program is successfully run with given inputs and outputs. If false, the next time running the program, skip_when_done=true will not take effect.\ndry_run::Bool = false: do not run the command, return (command::AbstractCmd, run_id_file::String).\n\nWorkflow\n\nValidate compatibility between p and inputs/outputs.\nCheck whether the program has run before. (skip_when_done, p.validate_outputs(outputs))\nCheck command dependencies. (check_dependencies, p.cmd_dependencies)\nValidate inputs. (p.validate_inputs(inputs))\nGenerate runnable command from p and inputs/outputs. (stdout, stderr, append)\nPreparing before running main command. (p.prerequisites(inputs, outputs))\nRun command generated in #5.\nValidate outputs. (p.validate_outputs(outputs))\nSuccess, touch run id file, and return true. (touch_run_id_file::Bool)\n\nExample\n\np = CmdProgram(\n\tcmd_dependencies = [julia],\n\tid_file = \"id_file\",\n\tinputs = [\"input\", \"input2\"],\n\toutputs = [\"output\"],\n\tcmd = `echo input input2 output`\n)\n\ninputs = Dict(\n\t\"input\" => `in1`,\n\t\"input2\" => 2\n)\n\noutputs = Dict(\n\t\"output\" => \"out\"\n)\n\nrun(p, inputs, outputs;\n\ttouch_run_id_file = false\n)\n\n\n\n\n\n","category":"method"},{"location":"API/#Pipelines.infer_outputs-Tuple{CmdProgram,Dict{String,V} where V}","page":"API","title":"Pipelines.infer_outputs","text":"infer_outputs(p::CmdProgram, inputs::Dict{String})\n\nInfer the default outputs from p::CmdProgram and inputs::Dict{String}.\n\n\n\n\n\n","category":"method"},{"location":"API/#Command-Dependency","page":"API","title":"Command Dependency","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"CmdDependency()\ncheck_dependency(::CmdDependency)\ncheck_dependency_dir(path::AbstractString; exit_when_false=true)\ncheck_dependency_file(path::AbstractString; exit_when_false=true)","category":"page"},{"location":"API/#Pipelines.CmdDependency-Tuple{}","page":"API","title":"Pipelines.CmdDependency","text":"Struct\n\nmutable struct CmdDependency\n\texec::Base.AbstractCmd\n\ttest_args::Base.AbstractCmd\n\tvalidate_success::Bool\n\tvalidate_stdout::Function\n\tvalidate_stderr::Function\n\texit_when_fail::Bool\nend\n\nMethods\n\nCmdDependency(;\n\texec::Base.AbstractCmd=``,\n\ttest_args::Base.AbstractCmd=``,\n\tvalidate_success::Bool=false,\n\tvalidate_stdout::Function=do_nothing,\n\tvalidate_stderr::Function=do_nothing,\n\texit_when_fail::Bool=true\n)\n\nCreate Command Dependency (CmdDependency).\n\nArguments\n\nexec::AbstractCmd: the command to call the dependency.\ntest_args::AbstractCmd: for testing purposes, the command to be appended to exec.\nvalidate_success::Bool: when checking the dependency, whether to validate the exit code == 0.\nvalidate_stdout::Function: a function takes standard out as String and return the validation result as ::Bool.\nvalidate_stderr::Function: a function takes standard error as String and return the validation result as ::Bool.\nexit_when_fail::Bool: if validation fails, whether to throw error and exit.\n\nExample\n\njulia = CmdDependency(\n\texec = Base.julia_cmd(),\n\ttest_args = `--version`,\n\tvalidate_success = true,\n\tvalidate_stdout = x -> occursin(r\"^julia version\", x),\n\tvalidate_stderr = do_nothing,\n\texit_when_fail = true\n)\n\ncheck_dependency(julia)\n\n\n\n\n\n","category":"method"},{"location":"API/#Pipelines.check_dependency-Tuple{CmdDependency}","page":"API","title":"Pipelines.check_dependency","text":"check_dependency(p::CmdDependency) -> Bool\n\nCheck CmdDependency by evaluating (pexec) (ptest_args).\n\nIf success, return true.\n\nIf fail, return false, or throw DependencyError when p.exit_when_fail set to true.\n\n\n\n\n\n","category":"method"},{"location":"API/#Pipelines.check_dependency_dir-Tuple{AbstractString}","page":"API","title":"Pipelines.check_dependency_dir","text":"check_dependency_dir(path::Union{AbstractString,Cmd}; exit_when_false=true) -> Bool\n\nChecke whether a directory exists. Return ::Bool.\n\n\n\n\n\n","category":"method"},{"location":"API/#Pipelines.check_dependency_file-Tuple{AbstractString}","page":"API","title":"Pipelines.check_dependency_file","text":"check_dependency_file(path::Union{AbstractString,Cmd}; exit_when_false=true) -> Bool\n\nChecke whether a file exists. Return ::Bool.\n\n\n\n\n\n","category":"method"},{"location":"API/#Utils","page":"API","title":"Utils","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"to_str(::Cmd)\nto_cmd(::Cmd)\nsplit(::Cmd)\nreplaceext(::String, ::AbstractString)\nremoveext(::String)","category":"page"},{"location":"API/#Pipelines.to_str-Tuple{Cmd}","page":"API","title":"Pipelines.to_str","text":"to_str(x) -> String\nstr(x) -> String\n\nConvert x to String.\n\nx::Cmd: remove backticks (return string(x)[2:end-1]).\nx::Nothing: return \"\".\nx::Vector: join elements with \"_\" as delim.\nx::Any: return string(x).\n\n\n\n\n\n","category":"method"},{"location":"API/#Pipelines.to_cmd-Tuple{Cmd}","page":"API","title":"Pipelines.to_cmd","text":"to_cmd(x) -> Cmd\n\nConvert x to Cmd.\n\nException: when x::Nothing, return nothing::Nothing.\n\n\n\n\n\n","category":"method"},{"location":"API/#Base.split-Tuple{Cmd}","page":"API","title":"Base.split","text":"split(c::Cmd) = c.exec :: Vector{String}\n\nReturn splitted arguments of Cmd.\n\n\n\n\n\n","category":"method"},{"location":"API/#Pipelines.replaceext-Tuple{String,AbstractString}","page":"API","title":"Pipelines.replaceext","text":"replaceext(path, replacement::AbstractString)\n\nIf the last component of a path contains a dot, leave everything before the dot as usual, and everything after the dot is replaced by replacement. Otherwise, replacement will be appended to path.\n\nIf replacement is empty, the last dot will be removed.\n\n\n\n\n\n","category":"method"},{"location":"API/#Pipelines.removeext-Tuple{String}","page":"API","title":"Pipelines.removeext","text":"removeext(path)\n\nIf the last component of a path contains a dot, leave everything before the dot as usual, and everything including and after the dot is discarded.\n\n\n\n\n\n","category":"method"}]
}
