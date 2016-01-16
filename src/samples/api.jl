### Measurement data
abstract MCData

### Model to estimate distribution from
abstract MCModel

### Information about the parameter to be estimated (default value and/or prior distribution)
abstract MCParameter

### Sampler types hold the components that fully specify a Monte Carlo sampler
abstract MCSampler

### Tuners of sampler hyper-parameters
abstract MCTuner

### Runner types indicate what type of simulation will be run (e.g. the Generalized MH way)
### Their fields fully specify the simulation details (e.g. total number or number of burn-in iterations)
abstract MCRunner

### Variable to hold the results of the simulation
abstract MCChain

### Types of the proposal distribution (e.g., Normal, Spline etc)
abstract ProposalDensity

###Derivative order of the sample
abstract DerivativeOrder
type NullOrder <: DerivativeOrder end
type FirstOrder <: DerivativeOrder end
type SecondOrder <: DerivativeOrder end
type ThirdOrder <: DerivativeOrder end
typealias GTNullOrder Union{FirstOrder,SecondOrder,ThirdOrder}
typealias GTFirstOrder Union{SecondOrder,ThirdOrder}

### Monte Carlo sample types hold a single Monte Carlo sample
### A typical Sample includes at least the vector of parameter values and the respective log-likelihood value
abstract MCSample{O<:DerivativeOrder}

### Heap types hold the temporary components used by a Monte Carlo sampler during its run
### This means that heap types represent the internal state ("local variables") of a Monte Carlo sampler
abstract MCSamplerState{S<:MCSample}

### Tune types hold the temporary output of the sampler that is used for tuning the sampler
abstract MCTunerState

typealias DefaultNumberType Float64

###Function to test that GMH module has loaded
function print_gmh_module_loaded()
  println("$module_name(current_module()) module loaded successfully")
end
