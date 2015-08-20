###ODE for the spring-mass dynamic system
###The second-order differentional equation
### a = -K/M * x
###is written as a system of first order differential equations
###Function arguments are as requested by the Sundials ODE solver package
### t = timepoints to evaluate the ODE
### y = the values of the variables (x,v=dx/dt)
### ydot = the derivate vaues (v=dx/dt,a=d2x/dt2)
### paras = the equation parameters (K,M)
@everywhere function spring_mass_ode(t,y,ydot,paras)
  ydot[1] = y[2]
  ydot[2] = -paras[1]/paras[2]*y[1] #-K/M*x
end

###Helper function to create a spring-mass MCModel
### initialconditions of the dynamic system
### defaultparams the default parameter values
### noisevar the measurement noise variance
### timepoints to evaluate the ODE, if not specified then a default range is used based on period of oscillation of the system
### priorfraction: parameter to determine the boundaries of the uniform priors of the parameters
function spring_mass_model(initialconditions,defaultparams,noisevar;timepoints =nothing,priorfraction =5.0)
  timepoints,measurements = spring_mass_data(initialconditions,defaultparams,noisevar;timepoints = timepoints)
  l,u = defaultparams-defaultparams./priorfraction,defaultparams+defaultparams./priorfraction #construct upper and lower bounds for the priors
  parameters = GeneralizedMetropolisHastings.ModelParameters(
    defaultparams,
    [Distributions.Uniform(l[1],u[1]),Distributions.Uniform(l[2],u[2])];
    names = ["K","M"])
  ODEModel(
    parameters,
    measurements,
    timepoints,
    repmat(noisevar,length(timepoints)),
    spring_mass_ode,
    initialconditions,2,[1,2];
    name = "Spring-Mass")
end
