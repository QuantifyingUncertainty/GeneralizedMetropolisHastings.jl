immutable ODEModel <: MCModel
  #Generic model specs
  name::AbstractString
  parameters::ModelParameters
  gradientepsilon::Float64

  #Data + noise
  measurements::Matrix{Float64}
  timepoints::Vector{Float64}
  variance::Matrix{Float64}

  #ODE specific
  ode::Function
  initial::Vector{Float64}
  nstates::Int
  observed::Vector{Int}
  unobserved::Vector{Int}
  abstol::Float64
  reltol::Float64
  inferinitial::Bool
end

###Constructor with default values for tolerance
ODEModel(parameters::ModelParameters,
         measurements::Matrix{Float64},
         timepoints::Vector{Float64},
         var::Matrix{Float64},
         ode::Function,
         initial::Vector{Float64},
         nstates::Int,
         observed::Vector{Int};
         name::AbstractString = "ODE",#From here on arguments are optional and have default values
         abstol::Float64 = 1e-6,
         reltol::Float64 = 1e-6,
         inferinitial = false) =
  ODEModel(name,parameters,100.0*abstol,measurements,timepoints,var,ode,initial,nstates,observed,setdiff(1:nstates,observed),abstol,reltol,inferinitial)

###Evaluate the model for the given parameters
evaluate(m::ODEModel,paras::Vector{Float64}) = Sundials.cvode((t,y,ydot)->m.ode(t,y,ydot,paras),m.initial,m.timepoints;reltol=m.reltol,abstol=m.abstol)[:,m.observed]

###Calculate the loglikelihood function / TODO: separate out the noise model
loglikelihood(m::ODEModel,sol::Matrix{Float64}) = sum(map((μ,σ,x)->logpdf(Distributions.Normal(μ,σ),x),m.measurements,sqrt(m.variance),sol))

###Calculate the gradient of the loglikelihood
gradienthelper(m::ODEModel,data::Matrix{Float64},sol::Matrix{Float64},grad::Array{Float64,3}) = vec(sum((grad.*(data-sol)./m.variance),(1,2))) #.* does automatic broadcast(), summing result along first 2 dimensions
gradloglikelihood(m::ODEModel,sol::Matrix{Float64},grad::Array{Float64,3}) = gradienthelper(m,m.measurements,sol,grad)

###Calculate the metric tensor
tensorvalue(m::ODEModel,grad::Array{Float64,3},i::Int,j::Int) = sum(grad[:,:,i].*grad[:,:,j]./m.variance)

###Generate pseudodata for the approximate metric tensor calculations
pseudodata(m::ODEModel,sol::Matrix{Float64}) = map((μ,σ)->rand(Distributions.Normal(μ,σ)),sol,sqrt(m.variance))
tangentvector(m::ODEModel,sol::Matrix{Float64},grad::Array{Float64,3}) = gradienthelper(m,pseudodata(m,sol),sol,grad)

###Base functionality
function Base.show(io::IO,m::ODEModel)
  println(io,"ODEModel ",m.name)
  show(io,m.parameters)
  println(io,"ODE function: ",m.ode)
  println(io,"number of measurements: ",size(m.measurements,1))
  println(io,"number of states: ",m.nstates)
  println(io,"observed: ",m.observed)
  println(io,"unobserved: ",m.unobserved)
  nothing
end
