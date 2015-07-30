@everywhere using GeneralizedMetropolisHastings
using Base.Test

###Helper function to calculate boundaries for a uniform prior around a certain parameter value
l_u(x::Float64,f::Float64) = x-x/f,x+x/f

###ODE for spring-mass dynamic system
@everywhere function spring_mass_ode(t,y,ydot,paras)
  ydot[1] = y[2]
  ydot[2] = -paras[1]/paras[2]*y[1] #-K/M*X
end

###Function that can generate date for the spring-mass dynamic system
function spring_mass_data(t,y0,paras,noisevar)
  a = sqrt(paras[1]/paras[2])
  y = [y0[1]*cos(a*t)+y0[2]/a*sin(a*t) -a*y0[1]*sin(a*t)+y0[2]*cos(a*t)]
  [y[:,1]+rand(Distributions.Normal(0.0,noisevar[1]),size(y,1)) y[:,2]+rand(Distributions.Normal(0.0,noisevar[2]),size(y,1))]
end

###Construct a model using the spring-mass ode and generated measurement data
function spring_mass_model(timepoints,initialconditions,defaultparams,variance)
  measurements = spring_mass_data(timepoints,initialconditions,defaultparams,variance)
  l1,u1 = l_u(defaultparams[1],5.0)
  l2,u2 = l_u(defaultparams[2],5.0)
  parameters::GeneralizedMetropolisHastings.ModelParameters = GeneralizedMetropolisHastings.ModelParameters(defaultparams,[Distributions.Uniform(l1,u1),Distributions.Uniform(l2,u2)]; names = ["K","M"])
  ODEModel(parameters,measurements,timepoints,repmat(variance,length(timepoints)),spring_mass_ode,initialconditions,2,[1,2]; name = "Spring + Mass")
end
