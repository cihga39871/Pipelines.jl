
mutable struct CmdDependency
    exec::Base.AbstractCmd
    test_args::Base.AbstractCmd
	validate_success::Bool
    validate_stdout::Function
    validate_stderr::Function
    exit_when_fail::Bool
end

function CmdDependency()
	CmdDependency(``, ``, false, do_nothing, do_nothing, true)
end

function CmdDependency(;
	exec::Base.AbstractCmd=``,
	test_args::Base.AbstractCmd=``,
	validate_success::Bool=false,
	validate_stdout::Function=do_nothing,
	validate_stderr::Function=do_nothing,
	exit_when_fail::Bool=true
	)
	CmdDependency(exec, test_args, validate_success, validate_stdout, validate_stderr, exit_when_fail)
end

"""
	readall(cmd::Cmd) -> (standard_out::String, standard_err::String, success::Bool)
"""
function readall(cmd::Base.AbstractCmd)
    out = Pipe()
    err = Pipe()

	process = try
	    run(pipeline(cmd, stdout=out, stderr=err))
	catch
		nothing
	end
    close(out.in)
    close(err.in)

    return (
        read(out, String),
        read(err, String),
        !isnothing(process)
    )
end

function check_dependency(p::CmdDependency)
    out, err, success = readall(`$(p.exec) $(p.test_args)`)
	has_dependency = true

	if p.validate_success && !success
		@goto dependency_error
	end

	try
		res = p.validate_stdout(out)
		res === false && error("DependencyError: invalid: $(p.exec): the `validate_stdout` function returns false.")
	catch e
		rethrow(e)
		@goto dependency_error
	end

	try
		res = p.validate_stderr(err)
		res === false && error("DependencyError: invalid: $(p.exec): the `validate_stdout` function returns false.")
	catch e
		rethrow(e)
		@goto dependency_error
	end

	return true

	@label dependency_error
	if p.exit_when_fail
		error("DependencyError: invalid: $(p.exec)")
	else
		@error "DependencyError: invalid: $(p.exec)"
	end
	return false
end

function Base.display(p::CmdDependency)
	print("CmdDependency\n  exec            :")
	display(p.exec)
	print("  test_args       :")
	display(p.test_args)
	println("  validate_success: $(p.validate_success)\n  validate_stdout : $(p.validate_stdout)\n  validate_stderr : $(p.validate_stderr)\n  exit_when_fail  : $(p.exit_when_fail)")
end

function Base.print(io::IO, p::CmdDependency)
	print(io, p.exec)
end

function Base.show(io::IO, p::CmdDependency)
	show(io, p.exec)
end


function check_dependency_file(path::AbstractString; exit_when_false=true)
	has_dependency = isfile(path)
	if exit_when_false && !has_dependency
		error("DependencyError: not a valid file: $path")
	end
	has_dependency
end
check_dependency_file(path::Cmd; exit_when_false=true) = check_dependency_file(str(path); exit_when_false=exit_when_false)

function check_dependency_dir(path::AbstractString; exit_when_false=true)
	has_dependency = isdir(path)
	if exit_when_false && !has_dependency
		error("DependencyError: not a valid directory: $path")
	end
	has_dependency
end
check_dependency_dir(path::Cmd; exit_when_false=true) = check_dependency_file(str(path); exit_when_false=exit_when_false)
