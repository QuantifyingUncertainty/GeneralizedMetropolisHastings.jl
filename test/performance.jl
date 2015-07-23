using GeneralizedMetropolisHastings
using Base.Test

performancetests = [
  "proposals",
  "samplers"]

println("==========================")
println("Running performance tests:")
println("==========================")

for t in performancetests
  tfile = t*".jl"
  println("  * $(tfile) *")
  include(string("performance/",tfile))
end


