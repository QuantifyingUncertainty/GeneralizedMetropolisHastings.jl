include("util.jl")
srand(0)

functionalitytests = [
  "springmasstest1"
  "springmasstest2"
  #"gaussiantest1",
  #"gaussiantest2"
  ]

println("============================")
println("Running functionality tests:")
println("============================")

for t in functionalitytests
  tfile = t*".jl"
  println("  * $(tfile) *")
  include(string("functionality/",tfile))
end

