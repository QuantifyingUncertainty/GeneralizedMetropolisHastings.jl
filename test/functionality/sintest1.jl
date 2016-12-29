###Test a Standard Metropolis-Hastings MCMC for a sine target model

###Innitialize variables

#MCMC iterations and burnin iterations
niterations = 100
nburnin = 100
ntunerperiod = 10

#Generalized MH with multiple proposals per iteration
nproposals = 30

#Time points to simulate the sine-cosine model
timepoints = linspace(0.0,10.0,100)

###Model parameter values (coefficients of sine and cosine)
alower = 2.0
areal = 3.0
aupper = 4.0

###The variance of the normal noise on the data
fcov = [0.01]

###Create a sine with measurement data and uniform prior on the parameter
fparas = parameters([:a],[alower],[aupper],[areal])
fdata = data(:array,timepoints,sin(3timepoints))
fnoise = noise(:gaussian,fcov)
fmodel = model(:target,fparas,fdata,fnoise,(p,t)->sin(p[1]*t),timepoints;name="SineTestModel")
show(fmodel)

###Create different samplers
fnparas = numparas(fmodel)
fsampler1 = sampler(:mh,:normal,0.1,eye(1))
fsampler2 = sampler(:adaptive,0.01,fnparas)

###Create different tuners for the samplers
ftuner1 = tuner(:scale,10,0.5,:erf)
ftuner2 = tuner(:monitor,10)

###Create a Generalized Metropolis-Hastings runner
p = policy(:mh,nproposals)
r = runner(p,niterations,nproposals;numburnin=nburnin)
show(r)

###Run the MCMC with MH sampler
c1 = run!(r,fmodel,fsampler1,ftuner1)
show(c1)

###Run the MCMC with Adaptive MH sampler
c2 = run!(r,fmodel,fsampler2,ftuner2)
show(c2)
