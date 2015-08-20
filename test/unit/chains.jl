import GeneralizedMetropolisHastings.MCChainOrdinary,
  GeneralizedMetropolisHastings.MCChainWithGradient,
  GeneralizedMetropolisHastings.create_chain,
  GeneralizedMetropolisHastings.store_results!,
  GeneralizedMetropolisHastings.store_common!,
  GeneralizedMetropolisHastings.store_gradient!,
  GeneralizedMetropolisHastings.add_acceptance!

s1 = MHNormal(2,1.0)
s2 = SmMALANormal(2,1.0)
h1 = MHHeap(s1,3)
h1.samples[1].values = rand(2)
h2 = SmMALAHeap(s2,3)
h2.samples[1].values = rand(2)

###Test the MCChain constructors
c1 = create_chain(2,10)
c2 = create_chain(2,10;storegradient = true,runtime = 1.0)
@test typeof(c1) == MCChainOrdinary
@test typeof(c2) == MCChainWithGradient
@test c1 == create_chain(2,10;storegradient = false)
@test c2.gradloglikelihood == zeros(2,10) && c2.runtime == 1.0

###Test add_acceptance
@test c1.acceptance == 0
add_acceptance!(c1,[1,1])
@test c1.acceptance == 0
add_acceptance!(c1,[1,2])
@test c1.acceptance == 1
add_acceptance!(c1,[2,1,3,3])
@test c1.acceptance == 3

###Test store_common!
c1.nsample = 1
@test c1.values[:,1] == zeros(2) && c1.loglikelihood[1] == 0.0 && c1.logprior[1] == 0.0
store_common!(c1,h1.samples[1])
@test c1.values[:,1] != zeros(2) && isequal(c1.loglikelihood[1],NaN) && isequal(c1.logprior[1],NaN)
c2.nsample = 1
@test c2.values[:,1] == zeros(2) && c2.loglikelihood[1] == 0.0 && c2.logprior[1] == 0.0 && c2.gradloglikelihood[:,1] == zeros(2)
store_common!(c2,h2.samples[1])
@test c2.values[:,1] != zeros(2) && isequal(c2.loglikelihood[1],NaN) && isequal(c2.logprior[1],NaN) && c2.gradloglikelihood[:,1] == zeros(2)
store_gradient!(c2,h2.samples[1])
@test isequal(c2.gradloglikelihood[:,1],fill(NaN,2))

###Test store_results!
c1 = create_chain(2,10)
c2 = create_chain(2,10;storegradient = true)
@test c1.nsample == 0
store_results!(c1,h1,[1,2,3])
@test c1.nsample == 2 && c1.acceptance == 2
@test c2.nsample == 0
store_results!(c2,h2,[1,2,2])
@test c2.nsample == 2 && c2.acceptance == 1

###Test show()
println()
println("====================")
println("Test show() function")
show(c1)
show(c2)
println("End  show() function")
println("====================")
println()
nothing
