###Test a Standard Metropolis-Hastings MCMC for a spring-mass ODE model

###Innitialize variables

#MCMC iterations and burnin iterations
niterations1 = 1000
nburnin1 = 500

#A single proposal makes it a Standard M-H (as opposed to Generalized with multiple proposals per iteration)
nproposals1 = 1

#Time points to simulate the spring-mass ODE
timepoints1 = 0:0.1:10

###Initial conditions for the spring-mass ODE (position and speed)
initial1 = [-1.0,1.0]

###Default values of the parameters (spring stiffness K and mass M)
K1 = 100.0 #in (N/m)
M1 = 10.0 #in kg

###The variance of the normal noise on the data
noise1 = [0.01,0.09]

###Create a spring-mass model with measurement data and ODE function
m1 = springmassmodel(timepoints1,initial1,[K1,M1],noise1,[K1-20.0,M1-2.0],[K1+20.0,M1+2.0])
show(m1)

###Create a Metropolis sampler with a Normal proposal density
s1 = sampler(:mh,:normal,0.1,[1.0 0.0;0.0 0.1])
show(s1)

###Create a tuner that scales the proposal density
t1 = tuner(:scale,100,0.5,:erf)
show(t1)

###Create a Generalized Metropolis-Hastings runner (which will default to Standard MH because nproposals == 1)
p1 = policy(:mh,nproposals1)
r1 = runner(p1,niterations1;numburnin=nburnin1)
show(r1)

###Run the MCMC
c1 = run!(r1,m1,s1,t1)
show(c1)
