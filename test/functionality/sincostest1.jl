###Test a Standard Metropolis-Hastings MCMC for a sin-cos target model

###Innitialize variables

#MCMC iterations and burnin iterations
niterations4 = 1000
nburnin4 = 500

#A single proposal makes it a Standard M-H (as opposed to Generalized with multiple proposals per iteration)
nproposals4 = 1

#Time points to simulate the sine-consine model
timepoints4 = linspace(0.0,10.0,200)

###Model parameter values (coefficients of sine and cosine)
a4 = 0.50
b4 = 2.0/3.0

###The variance of the normal noise on the data
noise4 = [0.01,0.04]

###Create a sine-cose model with measurement data and parameter priors
m4 = sincosmodel(timepoints4,[a4,b4],noise4,[0.1,0.1],[1.0,1.0])
show(m4)

###Create a Metropolis sampler with a Normal proposal density
s4 = sampler(:mh,:normal,0.001,eye(2))
show(s4)

###Create a tuner that scales the proposal density
t4 = tuner(:scale,100,0.5,:erf)
show(t4)

###Create a Generalized Metropolis-Hastings runner (which will default to Standard MH because nproposals == 1)
p4 = policy(:mh,nproposals4)
r4 = runner(p4,niterations4;numburnin=nburnin4)
show(r4)

###Run the MCMC
c4 = run!(r4,m4,s4,t4)
show(c4)
