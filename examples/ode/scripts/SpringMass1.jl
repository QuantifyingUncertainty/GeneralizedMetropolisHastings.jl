println("=================================================================")
println("Beginning execution of script examples/ode/scripts/SpringMass1.jl")
println("=================================================================")
###On a machine with more cores, you can increase NPROCS to use
###extra processes to speed up the Generalized Metropolis Hastings code
###This will initialize 1 main and NPROCS-1 worker processes
const NPROCS = 1
if nprocs() < NPROCS
    addprocs(NPROCS - nprocs())
end
println("Number of parallel processes running: ",nprocs())

###Add the src directory of the GeneralizedMetropolisHastings package to the path
push!(LOAD_PATH,"../../../src/")

###Make the GeneralizedMetropolisHastings code available to all processes
@everywhere using GeneralizedMetropolisHastings

###Use Pyplot to plot the results
using PyPlot

###Print a message indicating that the GMH package has loaded correctly
print_gmh_module_loaded()

println("================================")
println("Initialize Simulation Parameters")
println("================================")
#Total number of MCMC samples
nburnin = 1000
nsamples = 10000

#Standard M-H for nproposals == 1
#Generalized M-H for nproposals > 1
nproposals = 10

#Time points to simulate the spring-mass ODE
t = [0.0:0.1:10.0]

###Initial conditions for the spring-mass ODE (position and speed)
initialposition = -1.0 #in meters
initialvelocity = 1.0 #in meters/second

###Default values of the parameters (spring stiffness K and mass M)
K = 50.0 #in Newton/meter
M = 10.0 #in kg

###The variance of the normal noise on the input data
noise = [0.1 0.3]

println("==========================================")
println("Simulation parameters defined successfully")
println("==========================================")

###Make the data generation function available to the script
include("../data/springmass.jl")

###Make the ode and model creation function available to the script
include("../models/springmass.jl")

###Create a Spring-Mass model with measurement data and ODE function
m = spring_mass_model([initialposition,initialvelocity],[K,M],noise;timepoints =t)

###Plot the measurement data (simmulated data + noise)
figure("SpringMass1")
subplot(221)
plot(m.timepoints,m.measurements[:,1];label="location")
plot(m.timepoints,m.measurements[:,2];label="velocity")
xlabel("Time")
ylabel("Amplitude")
title("Spring-Mass measurement data")
grid("on")
legend(loc="upper right",fancybox="true")

###Show the model
println("==========================")
println("Model defined successfully")
println("==========================")
show(m)

###Create a Metropolis sampler with a Normal proposal density
s = MHNormal([5.0 0.0;0.0 1.0],0.1)
println("============================")
println("Sampler defined successfully")
println("============================")
show(s)

###Create a Generalized Metropolis-Hastings runner
###(which will default to Standard MH when nproposals == 1)
r = GMHRunner(nsamples,nproposals;burnin = nburnin)
println("===========================")
println("Runner defined successfully")
println("===========================")
show(r)

###Run the MCMC (can take quite a bit of time)
println("=======================")
println("Run the MCMC simulation")
println("=======================")
c = run!(r,m,s)
println("=========================")
println("Completed MCMC simulation")
println("=========================")
show(c)
println("Results of the MCMC simulation:")
println(" mean K: ",mean(c.values[1,nburnin:end]))
println(" mean M: ",mean(c.values[2,nburnin:end]))
println(" mean K/M: ",mean(c.values[1,nburnin:end]./c.values[2,nburnin:end]))
println("Mean K/M should be close to $(K/M)")

###Shut down the worker processes (to release resources before plotting)
for i in workers()
  if i != 1
    rmprocs(i)
  end
end

println("================")
println("Plotting results")
println("================")

###Plot the loglikelihood values across samples
###After an initial few low values, this should remain relatively high
subplot(222)
plot(1:length(c.loglikelihood),c.loglikelihood+c.logprior)
title("Loglikelihood values across samples")
xlabel("Samples")
ylabel("loglikelihood - logprior")

###Plot a scatter plot of K vs M values
###These should be spread around the K/M == 10.0 line (the diagonal in the figure)
ax3 = subplot(223)
kl,ku = K-K/5.0,K+K/5.0
ml,mu = M-M/5.0,M+M/5.0
ax3[:set_xlim]([kl,ku])
ax3[:set_ylim]([ml,mu])
scatter(c.values[1,:]',c.values[2,:]',marker=".",color="blue")
ax3[:set_aspect](abs(ku-kl)/abs(mu-ml))
title("MCMC samples of Spring-Mass ODE parameters")
xlabel("Stiffness K (N/m)")
ylabel("Mass M (kg)")
grid("on")

###Plot a histogram of K/M values, which should peak around the true ratio of K/M
ax4 = subplot(224)
kml,kmu = K/M-K/M/10.0,K/M+K/M/10.0
ax4[:set_xlim]([kml,kmu])
nbins = linspace(kml,kmu,100)
h = PyPlot.plt.hist((c.values[1,:]./c.values[2,:])',nbins)
grid("on")
xlabel("K/M")
ylabel("Number of Samples")
title("Histogram of K/M values")





