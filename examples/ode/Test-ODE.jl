###Use 8 processes (1 main and 7 workers) as JuliaBox gives you 8 processors
if nprocs() < 8
    addprocs(8-nprocs())
end
println("Number of parallel processes: ",nprocs())

push!(LOAD_PATH,"/home/juser/GeneralizedMetropolisHastings.jl/src/")

@everywhere using GeneralizedMetropolisHastings
using PyPlot

print_gmh_module_loaded()
