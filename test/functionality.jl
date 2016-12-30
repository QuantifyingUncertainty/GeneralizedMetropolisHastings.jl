testfolder = dirname(@__FILE__())

functionalitytests = [
    "sintest1"
  ]

println("============================")
println("Running functionality tests:")
println("============================")

for t in functionalitytests
  tfile = joinpath(testfolder,"functionality",string(t,".jl"))
  println("  * $(tfile) *")
  include(tfile)
  println()
end
