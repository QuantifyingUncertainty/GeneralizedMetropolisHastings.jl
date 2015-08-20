###Function to generate simulation data for the spring-mass dynamic system
### y0 = initial conditions [x,dx/dt]
### paras = [M,K]
### noisevar = variance of the measurement noise
### timepoints to evaluate the ODE; if not specified, use twice the period of oscillation
###Return: a tuple of:
### timepoints
### measurements
###  Column 1: Membrane Potential
###  Column 2: Recovery Variable
###Use examples/ode/test/plot_springmass.jl to plot data
function spring_mass_data(y0,paras,noisevar;timepoints =nothing)

  #the period of the spring-mass oscilation
  a = sqrt(paras[1]/paras[2])

  #generate timepoints (if not specified explicitly)
  if timepoints == nothing
    timepoints = linspace(0,4pi/a,100) #twice the period of oscillation
  end

  #generate the data
  y = [y0[1]*cos(a*timepoints)+y0[2]/a*sin(a*timepoints) -a*y0[1]*sin(a*timepoints)+y0[2]*cos(a*timepoints)]

  #add measurement noise
  y = [y[:,1]+rand(Distributions.Normal(0.0,noisevar[1]),size(y,1)) y[:,2]+rand(Distributions.Normal(0.0,noisevar[2]),size(y,1))]

  #return the data
  timepoints,y

end
