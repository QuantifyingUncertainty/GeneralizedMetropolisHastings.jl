### MCChain stores the output of a Monte Carlo iteration

type MCChain
  values::Array{Float64,3}
  loglikelihood::Matrix{Float64}
  gradloglikelihood::Array{Float64,3}
  diagnostics::Dict
  runtime::Float64
end

function MCChain(numparas::Int,numprops::Int,numiterations::Int;storegradient::Bool=false,diagnostics::Dict=Dict(),runtime::Float64=0.0)
  if storegradient
    MCChain(zeros(numparas,numprops,numiterations),zeros(numprops,numiterations),zeros(numparas,numprops,numiterations),diagnostics,runtime)
  else
    MCChain(zeros(numparas,numprops,numiterations),zeros(numprops,numiterations),zeros(0,0,0),diagnostics,runtime)
  end
end

function Base.show(io::IO,c::MCChain)
  nparas,nprops,niterations = size(c.values)
  println(io, "parameters: ",nparas,", proposals/iteration: ",niterations," $nsamples samples (per parameter), $(round(c.runtime, 1)) sec.")
end
