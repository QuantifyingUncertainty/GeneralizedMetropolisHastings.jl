###ODE for the FitzHugh-Nagumo model: http://www.scholarpedia.org/article/FitzHugh-Nagumo_model
### t = timepoints to evaluate the ODE at
### y = the values of the variables (Membrane Potential,Refractory Variable)
### ydot = the derivate vaues
### paras = the equation parameters (a,b,c)
@everywhere function fitzhugh_nagumo_ode(t,y,ydot,paras)

  a = paras[1]
  b = paras[2]
  c = paras[3]

  ydot[1] = c*(y[1]-(y[1]^3)/3+y[2])
  ydot[2] = -(y[1]-a+b*y[2])/c
end

###Helper function to create a FitzHugh-Nagumo MCModel
### initialconditions of the dynamic system
### defaultparams the default parameter values
### noisevar the measurement noise variance
function fitzhugh_nagumo_model(initialconditions,defaultparams,noisevar)
  timepoints,measurements = fitzhugh_nagumo_data()
  parameters = GeneralizedMetropolisHastings.ModelParameters(
    defaultparams,
    [Distributions.Uniform(0.0,5.0) for i=1:3];
    names = ["a","b","c"])
  ODEModel(parameters,
           measurements,
           timepoints,
           repmat(noisevar,length(timepoints)),
           fitzhugh_nagumo_ode,
           initialconditions,2,[1,2];
           name = "FitzHugh-Nagumo")
end
