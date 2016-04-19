include("imports.jl")
include("util.jl")
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
