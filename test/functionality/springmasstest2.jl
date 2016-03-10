###Test a Standard Metropolis-Hastings MCMC for a spring-mass ODE model

###Innitialize variables

#MCMC iterations and burnin iterations
niterations2 = 50
nburnin2 = 25

#Generalized M-H with multiple proposals per iteration
nproposals2 = 20

#Time points to simulate the spring-mass ODE
timepoints2 = 0:0.1:10

###Initial conditions for the spring-mass ODE (position and speed)
initial2 = [-1.0,1.0]

###Default values of the parameters (spring stiffness K and mass M)
K2 = 100.0 #in (N/m)
M2 = 10.0 #in kg

###The variance of the normal noise on the data
noise2 = [0.01,0.09]

###Create a spring-mass model with measurement data and ODE function
m2 = springmassmodel(timepoints2,initial2,[K2,M2],noise2,[K2-20.0,M2-2.0],[K2+20.0,M2+2.0])
show(m2)

###Create a Metropolis sampler with a Normal proposal density
s2 = sampler(:mh,:normal,0.1,[1.0 0.0;0.0 0.1])
show(s2)

###Create a tuner that scales the proposal density
t2 = tuner(:scale,5,0.5,:erf)
show(t2)

###Create a Generalized Metropolis-Hastings runner (which will run in as many samplerstates as there are processes running)
p2 = policy(:gmh,nproposals2)
r2 = runner(:gmh,niterations2,nproposals2,p2;numburnin=nburnin2)
show(r2)

###Run the MCMC
c2 = run!(r2,m2,s2;tuner=t2)
show(c2)
