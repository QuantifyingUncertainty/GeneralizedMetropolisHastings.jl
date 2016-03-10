srand(0) #to make the tests repeatable

unittests = [
    "policies",
    "parameters",
    "data",
    "samples",
    "proposals",
    "samplers",
    "indicators",
    "models",
    "tuners",
    "chains",
    "runners"]

println("===================")
println("Running unit tests:")
println("===================")

for t in unittests
    tfile = t*".jl"
    println("  * $(tfile) *")
    include(string("unit/",tfile))
    println()
    println()
end

