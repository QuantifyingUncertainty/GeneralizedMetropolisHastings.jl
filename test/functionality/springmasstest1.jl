###Test a Standard Metropolis-Hastings MCMC for a spring-mass ODE model

###Innitialize variables

#Total number of MCMC samples
nsamples = 2000

#A single proposal makes it a Standard M-H (as opposed to Generalized with multiple proposals per iteration)
nproposals = 1

#Time points to simulate the spring-mass ODE
timepoints = [0:0.1:10]

###Initial conditions for the spring-mass ODE (position and speed)
initial = [-1.0,1.0]

###Default values of the parameters (spring stiffness K and mass M)
K = 100.0 #in (N/m)
M = 10.0 #in kg

###The variance of the normal noise on the data
noise = [0.1 0.3]

###Create a spring-mass model with measurement data and ODE function
m = spring_mass_model(timepoints,initial,[K,M],noise)
show(m)

###Create a Metropolis sampler with a Normal proposal density
s = MHNormal([1.0 0.0;0.0 0.1],0.01)
show(s)

###Create a Generalized Metropolis-Hastings runner (which will default to Standard MH because nproposals == 1)
r = GMHRunner(nsamples,nproposals)
show(r)

###Run the MCMC
c = run!(r,m,s)
show(c)
