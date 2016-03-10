###Test a Standard Metropolis-Hastings MCMC for a spring-mass ODE model

###Innitialize variables

#MCMC iterations and burnin iterations
niterations3 = 50
nburnin3 = 25

#Generalized M-H with multiple proposals per iteration
nproposals3 = 24

#Time points to simulate the spring-mass ODE
timepoints3 = 0:0.1:10

###Initial conditions for the spring-mass ODE (position and speed)
initial3 = [-1.0,1.0]

###Default values of the parameters (spring stiffness K and mass M)
K3 = 100.0 #in (N/m)
M3 = 10.0 #in kg

###The variance of the normal noise on the data
noise3 = [0.01,0.09]

###Create a spring-mass model with measurement data and ODE function
m3 = springmassmodel(timepoints3,initial3,[K3,M3],noise3,[K3-20.0,M3-2.0],[K3+20.0,M3+2.0])
show(m3)

###Create a Metropolis sampler with a Normal proposal density
s3 = sampler(:mh,:normal,0.1,[1.0 0.0;0.0 0.1])
show(s3)

###Create a Generalized Metropolis-Hastings runner (which will run in 3 separate batches because of samplerstates=:test)
p3 = policy(:gmh,nproposals3;samplerstates=:test)
r3 = runner(:gmh,niterations3,nproposals3,p3;numburnin=nburnin3)
show(r3)

###Run the MCMC
c3 = run!(r3,m3,s3)
show(c3)
