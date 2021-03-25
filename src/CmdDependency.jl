
mutable struct CmdDependency
    exec::Base.AbstractCmd
    test_args::Base.AbstractCmd
	validate_success::Bool
    validate_stdout::Function
    validate_stderr::Function
    exit_when_fail::Bool
end

"""
# Struct

	mutable struct CmdDependency
		exec::Base.AbstractCmd
		test_args::Base.AbstractCmd
		validate_success::Bool
		validate_stdout::Function
		validate_stderr::Function
		exit_when_fail::Bool
	end

# Methods

	CmdDependency(;
		exec::Base.AbstractCmd=``,
		test_args::Base.AbstractCmd=``,
		validate_success::Bool=false,
		validate_stdout::Function=do_nothing,
		validate_stderr::Function=do_nothing,
		exit_when_fail::Bool=true
	)

Create Command Dependency (`CmdDependency`).

# Arguments

- `exec::AbstractCmd`: the command to call the dependency.

- `test_args::AbstractCmd`: for testing purposes, the command to be appended to `exec`.

- `validate_success::Bool`: when checking the dependency, whether to validate the exit code == 0.

- `validate_stdout::Function`: a function takes standard out as `String` and return the validation result as `::Bool`.

- `validate_stderr::Function`: a function takes standard error as `String` and return the validation result as `::Bool`.

- `exit_when_fail::Bool`: if validation fails, whether to throw error and exit.

# Example

	julia = CmdDependency(
		exec = Base.julia_cmd(),
		test_args = `--version`,
		validate_success = true,
		validate_stdout = x -> occursin(r"^julia version", x),
		validate_stderr = do_nothing,
		exit_when_fail = true
	)

	check_dependency(julia)

"""
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

"""
	check_dependency(p::CmdDependency) -> Bool

Check `CmdDependency` by evaluating `\`\$(p.exec) \$(p.test_args)\``.

If success, return `true`.

If fail, return `false`, or throw DependencyError when `p.exit_when_fail` set to `true`.
"""
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

"""
	check_dependency_file(path::Union{AbstractString,Cmd}; exit_when_false=true) -> Bool

Checke whether a file exists. Return `::Bool`.
"""
function check_dependency_file(path::AbstractString; exit_when_false=true)
	has_dependency = isfile(path)
	if exit_when_false && !has_dependency
		error("DependencyError: not a valid file: $path")
	end
	has_dependency
end
check_dependency_file(path::Cmd; exit_when_false=true) = check_dependency_file(str(path); exit_when_false=exit_when_false)


"""
	check_dependency_dir(path::Union{AbstractString,Cmd}; exit_when_false=true) -> Bool

Checke whether a directory exists. Return `::Bool`.
"""
function check_dependency_dir(path::AbstractString; exit_when_false=true)
	has_dependency = isdir(path)
	if exit_when_false && !has_dependency
		error("DependencyError: not a valid directory: $path")
	end
	has_dependency
end
check_dependency_dir(path::Cmd; exit_when_false=true) = check_dependency_file(str(path); exit_when_false=exit_when_false)

"""
	exec(p::CmdDependency) -> AbstractCmd

Return the command to call the dependency.
"""
exec(p::CmdDependency) = p.exec