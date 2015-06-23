using GeneralizedMetropolisHastings
using Base.Test

tests =
  ["odetest1",
   "gaussiantest1",
   "gaussiantest2"]

println("Running tests:")

for t in tests
  tfile = t*".jl"
  println("  * $(tfile) *")
  include(tfile)
end
