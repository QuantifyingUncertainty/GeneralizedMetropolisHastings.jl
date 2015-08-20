@everywhere using GeneralizedMetropolisHastings
using Base.Test

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

###Helper function to create a spring-mass MCModel
### timepoints to evaluate the ODE
### initialconditions of the dynamic system
### defaultparams the default parameter values
### noisevar the measurement noise variance
### priorfraction: parameter to determine the boundaries of the uniform priors of the parameters
function spring_mass_model(timepoints,initialconditions,defaultparams,noisevar,priorfraction =5.0)
  measurements = spring_mass_data(timepoints,initialconditions,defaultparams,noisevar)
  l,u = defaultparams-defaultparams./priorfraction,defaultparams+defaultparams./priorfraction #construct upper and lower bounds for the priors
  parameters = GeneralizedMetropolisHastings.ModelParameters(defaultparams,[Distributions.Uniform(l[1],u[1]),Distributions.Uniform(l[2],u[2])]; names = ["K","M"])
  ODEModel(parameters,measurements,timepoints,repmat(noisevar,length(timepoints)),spring_mass_ode,initialconditions,2,[1,2]; name = "Spring-Mass")
end
