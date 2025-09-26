using Test
using Pipelines

@testset "redirect" begin

    @show stdout
    @show stderr
    @test stdout isa ScopedStreams.ScopedStream
    @test stderr isa ScopedStreams.ScopedStream
    
    ## cmd program
    p = CmdProgram(
        id_file = "asv",
        inputs = [
            :a => 10.6 => Float64,
            :b =>  5 => Int
        ],
        outputs = "c" => "<a>.<b>",
        cmd = `echo this is stdout a b` & pipeline(`julia -e '@info "this is stderr"'`)
    )
    out_io = open("out.txt", "w+")

    @test run(p, touch_run_id_file=false, stdout=out_io, stderr="err.txt", stdlog="log.txt")[1]

    close(out_io)

    @test read("out.txt", String) == "this is stdout 10.6 5\n"
    @test read("err.txt", String) == "[ Info: this is stderr\n"
    @test !occursin("[ Info: this is stderr", read("log.txt", String))

    rm("out.txt", force=true)
    rm("err.txt", force=true)
    rm("log.txt", force=true)

    p2 = CmdProgram(
        id_file = "vsad",
        inputs = [
            :a => 10.6 => Float64,
            :b =>  5 => Int
        ],
        outputs = "c" => "<a>.<b>",
        cmd = pipeline(`echo 123`, "out2.txt")
    )
    @test run(p2, touch_run_id_file=false, stdout=stdout)[1]  # stdout=stdout will not take effect because redirection was done in out2.txt
    @test isfile("out2.txt") && read("out2.txt", String) == "123\n"
    rm("out2.txt", force=true)
end