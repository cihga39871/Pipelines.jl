
### stdout and stderr redirecting

stdout_origin = nothing  # re-defined in __init__()
stderr_origin = nothing  # re-defined in __init__()

"""
    restore_stdout()

Restore the current stdout to the original stdout. It is useful when redirecting stdout/stderr fails when calling `redirect_to_files`, which happens when an old stream is closed and then redirected to.

!!! warning "Thread safety"
    Redirecting in Julia are not thread safe, so unexpected redirection might be happen if you are running programs in different `Tasks` or multi-thread mode.

See also `restore_stderr()`.
"""
function restore_stdout()
    global stdout_origin
    if isnothing(stdout_origin)
        @error "Failed to redirect since original stdout is not stored."
    else
        redirect_stdout(stdout_origin)
    end
end

"""
    restore_stderr()

Restore the current stderr to the original stderr. It is useful when redirecting stdout/stderr fails when calling `redirect_to_files`, which happens when an old stream is closed and then redirected to.

!!! warning "Thread safety"
    Redirecting in Julia are not thread safe, so unexpected redirection might be happen if you are running programs in different `Tasks` or multi-thread mode.

See also `restore_stdout()`.
"""
function restore_stderr()
    global stderr_origin
    if isnothing(stderr_origin)
        @error "Failed to redirect since original stderr is not stored."
    else
        redirect_stderr(stderr_origin)
    end
end



### capture and display StackTraceVector

"""
    struct StackTraceVector
        x::Vector
    end

- `x = [(exception,backtrace), ...]`: the result of `Base.current_exceptions()` in Julia 1.7 or `Base.catch_stack()` in Julia 1.1-1.6.
"""
struct StackTraceVector
    x::Vector
end
function Base.show(io::IO, v::StackTraceVector)
    for (exc, bt) in v.x
        showerror(io, exc, bt)
        println(io)
    end
end
function Base.show(v::StackTraceVector)
    for (exc, bt) in v.x
        showerror(stderr, exc, bt)
        println(stderr)
    end
end
function Base.display(v::StackTraceVector)
    for (exc, bt) in v.x
        showerror(stderr, exc, bt)
        println(stderr)
    end
end


"""
    try_function(f::Function, error_io::IO)
    try_function(f::Function, ::Nothing   )

Try to run `f`. If `f` throws error, display stacktraces in `error_io` or `stderr`, and return stacktrace information stored as `::StackTraceVector`.
"""
function try_function end

@static if VERSION < v"1.1"
    # no catch_stack
    function try_function(f::Function, io::IO)
        logger = SimpleLogger(io)
        res = with_logger(logger) do
            try
                f()
            catch exc
                bt = catch_backtrace()
                stv = StackTraceVector([(exc, bt)])
                show(io, stv)
                stv
            end
        end
    end
    function try_function(f::Function, ::Nothing)
        res = try
            f()
        catch exc
            bt = catch_backtrace()
            stv = StackTraceVector([(exc, bt)])
            show(stv)
            stv
        end
    end
elseif VERSION < v"1.7"
    # function named as catch_stack
    function try_function(f::Function, io::IO)
        logger = SimpleLogger(io)
        res = with_logger(logger) do
            try
                f()
            catch
                exps = StackTraceVector(Base.catch_stack())
                show(io, exps)
                exps
            end
        end
    end
    function try_function(f::Function, ::Nothing)
        res = try
            f()
        catch
            exps = StackTraceVector(Base.catch_stack())
            show(exps)
            exps
        end
    end
else
    # function named as current_exceptions after 1.7
    function try_function(f::Function, io::IO)
        logger = SimpleLogger(io)
        res = with_logger(logger) do
            try
                f()
            catch
                StackTraceVector(Base.current_exceptions()) # show(exps)
                # show(io, exps)
                # exps
            end
        end
    end
    function try_function(f::Function, ::Nothing)
        res = try
            f()
        catch
            StackTraceVector(Base.current_exceptions()) # show(exps)
            # show(exps)
            # exps
        end
    end
end



### extend IO

Base.close(::Nothing) = nothing
# Base.SimpleLogger(::Nothing) = nothing
# Base.with_logger(f::Function, ::Nothing) = f()

Base.redirect_stdout(f::Function, ::Nothing) = f()
# redirect_stdout and redirect_stderr are the same (::Base.RedirectStdStream) at least from julia v1.7, so defining redirect_stdout means redirect_stderr is also defined. If diff exists in previous julia versions, check first.
try
    Base.redirect_stderr(time, nothing)
catch
    Base.redirect_stderr(f::Function, ::Nothing) = f()
end


handle_open(::Nothing, mode) = nothing
handle_open(io::IO, mode) = io # do not change and do not close when exit
handle_open(file::AbstractString, mode) = open(file::AbstractString, mode)



### main redirection

"""
    redirect_to_files(f::Function, file; mode="a+")
    redirect_to_files(f::Function, outfile, errfile; mode="a+")
    redirect_to_files(f::Function, outfile, errfile, logfile; mode="a+")

Redirect outputs of function `f` to file(s).

- `xxxfile`: File path (`AbstractString`), `nothing` or `::IO`. `nothing` means no redirect. Files can be the same.
- `mode`: same as `open(..., mode)`.

Caution: If `xxxfile` is an `IO`, it won't be closed. Please use `close(io)` or `JobSchedulers.close_in_future(io, jobs)` manually!

!!! warning "Thread safety"
    Redirecting in Julia are not thread safe, so unexpected redirection might be happen if you are running programs in different `Tasks` or multi-thread mode.
"""
function redirect_to_files(f::Function, outfile, errfile, logfile; mode="a+")
    out = handle_open(outfile, mode)
    err = errfile == outfile ? out : handle_open(errfile, mode)
    log = logfile == outfile ? out : logfile == errfile ? err : handle_open(logfile, mode)

    old_stdout = Base.stdout
    old_stderr = Base.stderr

    if !isnothing(out)
        try
            redirect_stdout(out)
        catch
        end
    end
    if !isnothing(err)
        try
            redirect_stderr(err)
        catch
        end
    end

    # show_error = Base.stderr !== stderr_origin

    res = try_function(f, log)

    # before switch to old stdxxx, check whether it is opened. It might be happen due to redirect_xxx is not thread safe!
    if !isnothing(out)
        if isopen(old_stdout)
            try
                # it is not atomic, so use try
                redirect_stdout(old_stdout)
            catch
                try
                    redirect_stdout(stdout_origin)
                catch
                end
            end
        else
            try
                redirect_stdout(stdout_origin)
            catch
            end
        end
    end

    if !isnothing(err)
        if isopen(old_stderr)
            try
                # it is not atomic, so use try
                redirect_stderr(old_stderr)
            catch
                try
                    redirect_stderr(stderr_origin)
                catch
                end
            end
        else
            try
                redirect_stderr(stderr_origin)
            catch
            end
        end
    end

    outfile isa IO ? flush(outfile) : close(out)
    errfile isa IO ? flush(errfile) : close(err)
    logfile isa IO ? flush(logfile) : close(log)
    # if res isa StackTraceVector && show_error
    #     show(res)
    # end
    return res
end

redirect_to_files(f::Function, outfile, errfile; mode="a+") = redirect_to_files(f::Function, outfile, errfile, errfile; mode=mode)

redirect_to_files(f::Function, outfile; mode="a+") = redirect_to_files(f::Function, outfile, outfile, outfile; mode=mode)
