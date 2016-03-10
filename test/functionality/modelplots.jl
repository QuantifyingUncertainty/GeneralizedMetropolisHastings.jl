using PyPlot

#Time points to simulate the spring-mass ODE
timepoints1 = 0:0.1:10

###Initial conditions for the spring-mass ODE (position and speed)
initial = [-1.0,1.0]

###Default values of the parameters (spring stiffness K and mass M)
K = 100.0 #in (N/m)
M = 10.0 #in kg

###The variance of the normal noise on the data
noise1 = [0.01,0.09]

###Create a spring-mass model with measurement data and ODE function
m1 = springmassmodel(timepoints1,initial,[K,M],noise1,[K,M])
show(m1)

######################

#Time points to simulate the sine-cosine model
timepoints2 = linspace(0.0,10.0,200)

###The variance of the normal noise on the data
noise2 = [0.01,0.04]

###Create the sine-cosine model
m2 = sincosmodel(timepoints2,[0.5,2/3],noise2,[1.0,1.0])
show(m2)

######################

figure("functionality/modelplots")
subplot(211)
plot(timepoints1,datavalues(m1.measurements))
title(m1.name)
subplot(212)
title(m2.name)
plot(timepoints2,datavalues(m2.measurements))
