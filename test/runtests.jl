using GeneralizedMetropolisHastings
using Base.Test

unittests = [
  "indicatorsampletest"]

functionalitytests = [
  "odetest1",
  "gaussiantest1",
  "gaussiantest2"]

println("==================")
println("Running unittests:")
println("==================")

for t in unittests
  tfile = t*".jl"
  println("  * $(tfile) *")
  include(tfile)
end

println("===========================")
println("Running functionalitytests:")
println("===========================")

# for t in functionalitytests
#   tfile = t*".jl"
#   println("  * $(tfile) *")
#   include(tfile)
# end
