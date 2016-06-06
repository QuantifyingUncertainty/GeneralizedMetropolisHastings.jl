srand(0)

functionalitytests = [
    "distributions",
    "modelplots",
    "springmasstest1",
    "springmasstest2",
    "springmasstest3",
    "sincostest1",
    "sincostest2"
  ]

println("============================")
println("Running functionality tests:")
println("============================")

for t in functionalitytests
  tfile = t*".jl"
  println("  * $(tfile) *")
  include(string("functionality/",tfile))
end

