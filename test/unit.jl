using GeneralizedMetropolisHastings
using Base.Test

unittests = [
  #"indicatorsampletest",
  "parameters",
  "policies",
  "samples",
  "proposals",
  "samplers"]

println("===================")
println("Running unit tests:")
println("===================")

for t in unittests
  tfile = t*".jl"
  println("  * $(tfile) *")
  include(string("unit/",tfile))
end

