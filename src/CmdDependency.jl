
mutable struct CmdDependency
    exec::Base.Cmd
    test_args::Base.Cmd
    validate_success::Bool
    validate_stdout::Function
    validate_stderr::Function
    exit_when_fail::Bool
    status::UInt8
end

const NOT_CHECKED = 0b00
const OK = 0b01
const FAILED = 0b10

"""
# Struct

    mutable struct CmdDependency
        exec::Base.Cmd
        test_args::Base.Cmd
        validate_success::Bool
        validate_stdout::Function
        validate_stderr::Function
        exit_when_fail::Bool
        status::UInt8
    end

# Methods

    CmdDependency(;
        exec::Base.Cmd=``,
        test_args::Base.Cmd=``,
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
    exec::Cmd=``,
    test_args::Cmd=``,
    validate_success::Bool=false,
    validate_stdout::Function=do_nothing,
    validate_stderr::Function=do_nothing,
    exit_when_fail::Bool=true
    )
    CmdDependency(exec, test_args, validate_success, validate_stdout, validate_stderr, exit_when_fail, NOT_CHECKED)
end

"""
    readall(cmd::Cmd) -> (standard_out::String, standard_err::String, success::Bool)
"""
function readall(cmd::Base.AbstractCmd)
    out = Pipe()
    err = Pipe()

    process = try
        run(cmd, devnull, out, err)
    catch
        nothing
    end
    close(out.in)
    close(err.in)

    standard_out = read(out, String)
    standard_err = read(err, String)
    close(out.out)
    close(err.out)

    return (
        standard_out,
        standard_err,
        !isnothing(process)
    )
end

"""
    check_dependency(p::CmdDependency; exit_when_fail::Bool = p.exit_when_fail) -> Bool

Check `CmdDependency` by evaluating:

    `\$(p.exec) \$(p.test_args)`

If success, return `true`.

If fail, return `false`, or throw DependencyError when `exit_when_fail` set to `true`.
"""
function check_dependency(p::CmdDependency; exit_when_fail::Bool = p.exit_when_fail)

    if p.status == OK
        return true
    elseif p.status == FAILED
        @goto quick_fail # COV_EXCL_LINE
    end

    out, err, success = readall(`$p $(p.test_args)`)

    if p.validate_success && !success
        p.status = FAILED
        @label quick_fail # COV_EXCL_LINE
        if exit_when_fail
            error("CmdDependencyError: invalid: $p with $(p.test_args)")
        else
            @error timestamp() * "CmdDependencyError: invalid: $p $(p.test_args)" _module=nothing _group=nothing _id=nothing _file=nothing
        end
        return false
    end

    res1 = isok(p.validate_stdout(out))
    res2 = isok(p.validate_stderr(err))

    if res1 == res2 == true
        p.status = OK
        return true
    end

    p.status = FAILED
    error_msg = replace("CmdDependencyError: invalid: $p with $(p.test_args)\n   stdout: $out\n   stderr: $err", r"\n+" => "\n")
    if exit_when_fail
        error(error_msg)
    else
        @error timestamp() * error_msg CHECK_ARGS=p.test_args _module=nothing _group=nothing _id=nothing _file=nothing
    end
    return false
end


# Interpolation in Cmd
# It allows CmdDependency to be interpolated in `$p`.
Base.arg_gen(p::CmdDependency) = Base.arg_gen(p.exec)

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
