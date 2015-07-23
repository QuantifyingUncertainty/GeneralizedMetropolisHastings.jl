using GeneralizedMetropolisHastings
using Base.Test

alltests = [
  #"indicatorsampletest",
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
