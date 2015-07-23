using GeneralizedMetropolisHastings
using Base.Test

functionalitytests = [
  "odetest1",
  "gaussiantest1",
  "gaussiantest2"]

println("============================")
println("Running functionality tests:")
println("============================")

for t in functionalitytests
  tfile = t*".jl"
  println("  * $(tfile) *")
  include(string("functionality/",tfile))
end

