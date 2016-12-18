srand(0) #to make the tests repeatable

unittests = [
    "policies",
    "parameters",
    "distributions",
    "data",
    "samples",
    "proposals",
    "samplers",
    "indicators",
    "models",
    "tuners",
    "chains",
    "jobsegments",
    "runners",
    "mhrunner",
    "smhrunner",
    "gmhrunner"]

println("===================")
println("Running unit tests:")
println("===================")

include("unit/testutil.jl")

for t in unittests
    tfile = t*".jl"
    println("  * $(tfile) *")
    include(string("unit/",tfile))
    println()
    println()
end


