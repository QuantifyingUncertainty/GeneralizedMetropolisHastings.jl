println("=====================================================================")
println("Beginning execution of script examples/ode/scripts/FitzHughNagumo1.jl")
println("=====================================================================")
###Use extra processes to run the Generalized Metropolis Hastings code
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

###Initial conditions for the spring-mass ODE (membrane potential and refractory variable)
initialconditions = [-1.0,1.0]

###Default values of the parameters (a,b,c)
defaultparams = [0.2,0.2,3.0]

###The variance of the noise on the input data
noise = [0.02 0.005]

println("==========================================")
println("Simulation parameters defined successfully")
println("==========================================")

###Make the data generation function available to the script
include("../data/fitzhughnagumo.jl")

###Make the ode and model creation function available to the script
include("../models/fitzhughnagumo.jl")

###Create a FitzHugh-Nagumo model with measurement data and ODE function
m = fitzhugh_nagumo_model(initialconditions,defaultparams,noise)

###Plot the measurement data (simmulated data + noise)
figure("FitzHughNagumo1")
subplot(231)
plot(m.timepoints,m.measurements[:,1];label="membrane potential")
plot(m.timepoints,m.measurements[:,2];label="refractory variable")
xlabel("Time")
ylabel("Amplitude")
title("FitzHugh-Nagumo measurement data")
grid("on")
legend(loc="upper right",fancybox="true")

###Show the model
println("==========================")
println("Model defined successfully")
println("==========================")
show(m)

###Create a Metropolis sampler with a Normal proposal density
s = MHNormal(eye(3),0.01)
println("============================")
println("Sampler defined successfully")
println("============================")
show(s)

###Create a Generalized Metropolis-Hastings runner (which will default to Standard MH because nproposals == 1)
r = GMHRunner(nsamples,nproposals;burnin = nburnin,initialize =ValuesFromDefault())
println("===========================")
println("Runner defined successfully")
println("===========================")
show(r)

###Run the MCMC (can take quite a bit of time)
println("=======================")
println("Run the MCMC simulation")
println("=======================")
c = run!(r,m,s)
println("==========================")
println("Completed MCMC simulation")
println("=========================")

###Shut down the worker processes (to release resources before plotting)
for i in workers()
  if i != 1
    rmprocs(i)
  end
end

###Show the results of the simulations
show(c)
println("Results of the MCMC simulation:")
for i=1:3
  println(" parameter $(m.parameters.names[i]):  mean = ",mean(c.values[i,nburnin:end])," std = ",std(c.values[i,nburnin:end]))
end

println("================")
println("Plotting results")
println("================")

###Plot the loglikelihood values across samples
###After an initial few low values, this should remain relatively high
subplot(232)
plot(1:length(c.loglikelihood),c.loglikelihood+c.logprior)
title("Loglikelihood values across samples")
xlabel("Samples")
ylabel("loglikelihood - logprior")

###Plot the histograms of a,b,c values
for i=1:3
  subplot(233 + i)
  h = PyPlot.plt.hist(c.values[i,:]',20)
  grid("on")
  title("Parameter $(m.parameters.names[i])")
  xlabel("Values")
  ylabel("Samples")
end







