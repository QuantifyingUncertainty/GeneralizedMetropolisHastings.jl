srand(0)

performancetests = []

println("==========================")
println("Running performance tests:")
println("==========================")

for t in performancetests
  tfile = t*".jl"
  println("  * $(tfile) *")
  include(string("performance/",tfile))
end


