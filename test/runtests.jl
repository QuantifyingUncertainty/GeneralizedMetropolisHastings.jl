#run everything on one process
if nprocs() > 1
    rmprocs(workers())
end

println()
println("===========================================")
println("+++++++++++++++++++++++++++++++++++++++++++")
println("Running all tests with ",nprocs()," process")
println("+++++++++++++++++++++++++++++++++++++++++++")
println("===========================================")
println()

#import all functionality
include("imports.jl")
srand(0)

#all test categories
alltests = [
    "unit",
    "performance",
    "functionality"]

println("==================")
println("Running all tests:")
println("==================")

for t in alltests
    tfile = t*".jl"
    println("  * $(tfile) *")
    include(tfile)
end

#run the runner tests on additional processes
addprocs(2)
println()
println("================================================")
println("++++++++++++++++++++++++++++++++++++++++++++++++")
println("Running the multiprocess code tests with ",nprocs()," processes")
println("++++++++++++++++++++++++++++++++++++++++++++++++")
println("================================================")
println()

for i=1:nworkers()
    remotecall_wait(include,workers()[i],"test/unit/testutil.jl")
end

partests = [
    "unit/jobsegments",
    "unit/runners",
    "unit/mhrunner",
    "unit/smhrunner",
    "unit/gmhrunner",
    "unit/multiprocess",
    "functionality/sintest1"]

for t in partests
    tfile = t*".jl"
    println("  * $(tfile) *")
    try
        include(tfile)
    catch e
        if nprocs() > 1
            rmprocs(workers())
        end
        throw(e)
    end
end

if nprocs() > 1
    rmprocs(workers())
end
sleep(1.0)
