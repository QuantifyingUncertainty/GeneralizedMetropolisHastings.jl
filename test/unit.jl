testfolder = dirname(@__FILE__())

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

for t in unittests
    tfile = joinpath(testfolder,"unit",string(t,".jl"))
    println("  * $(tfile) *")
    include(tfile)
    println()
end
