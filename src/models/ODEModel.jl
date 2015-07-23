
type ODEModel <: MCModel

  #generic model
  name::String
  parameters::ModelParameters

  # Data
  measurements::Matrix{Float64}
  timepoints::Vector{Float64}
  variance::Matrix{Float64}

  # ODE specific
  ode::Function
  initial::Vector{Float64}
  nstates::Int
  observed::Vector{Int}
  unobserved::Vector{Int}
  abstol::Float64
  reltol::Float64
  inferinitial::Bool

end

evaluate(m::ODEModel,paras::Vector{Float64}) = Sundials.cvode((t,y,ydot)->m.ode(t,y,ydot,paras),m.initial,m.timepoints[:,1];reltol=m.reltol,abstol=m.abstol)[:,m.observed]
loglikelihood(m::ODEModel,sol::Matrix{Float64}) = sum(map((μ,σ,x)->logpdf(Normal(μ,σ),x),m.measurements,sqrt(m.variance),sol))
