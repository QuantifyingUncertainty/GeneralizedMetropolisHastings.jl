###Test a Standard Metropolis-Hastings MCMC for a sin-cos target model

###Innitialize variables

#MCMC iterations and burnin iterations
niterations5 = 100
nburnin5 = 10

#Generalized MH with multiple proposals per iteration
nproposals5 = 50

#Time points to simulate the sine-cosine model
timepoints5 = linspace(0.0,10.0,200)

###Model parameter values (coefficients of sine and cosine)
a5 = 0.50
b5 = 2.0/3.0

###The variance of the normal noise on the data
noise5 = [0.01,0.04]

###Create a sine-cosine with measurement data and uniform priors on the parameters
m5 = sincosmodel!(timepoints5,[a5,b5],noise5,[0.1,0.1],[1.0,1.0])
show(m5)

###Create a Metropolis sampler with a Normal proposal density
s5 = sampler(:mh,:normal,0.001,eye(2))
show(s5)

###Create a tuner that scales the proposal density
t5 = tuner(:scale,2,0.5,:erf)
show(t5)

###Create a Generalized Metropolis-Hastings runner (which will default to Standard MH because nproposals == 1)
p5 = policy(:gmh,nproposals5)
r5 = runner(:gmh,niterations5,nproposals5,p5;numburnin=nburnin5)
show(r5)

###Run the MCMC
c5 = run!(r5,m5,s5;tuner=t5)
show(c5)
