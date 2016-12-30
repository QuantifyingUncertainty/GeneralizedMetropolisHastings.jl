#run everything on one process
if nprocs() > 1
    rmprocs(workers())
end

#get the correct paths and filenames on all systems
testfolder = dirname(@__FILE__())
importsfile = joinpath(testfolder,"imports.jl")
utilfile = joinpath(testfolder,"testutil.jl")

#import all functionality on the local process
include(importsfile)
include(utilfile)

#make the run repeatable
srand(0)

println()
println("===========================================")
println("+++++++++++++++++++++++++++++++++++++++++++")
println("Running all tests with ",nprocs()," process")
println("+++++++++++++++++++++++++++++++++++++++++++")
println("===========================================")
println()

#all test categories
alltests = [
    "unit",
    "performance",
    "functionality"]

println("==================")
println("Running all tests:")
println("==================")

for t in alltests
    tfile = joinpath(testfolder,string(t,".jl"))
    println("  * $(tfile) *")
    include(tfile)
    println()
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
    remotecall_wait(include,workers()[i],importsfile)
    remotecall_wait(include,workers()[i],utilfile)
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
    tfile = joinpath(testfolder,string(t,".jl"))
    println("  * $(tfile) *")
    try
        include(tfile)
    catch e
        if nprocs() > 1
            rmprocs(workers())
        end
        throw(e)
    end
    println()
end

if nprocs() > 1
    rmprocs(workers())
end
sleep(1.0)
