@test GeneralizedMetropolisHastings._numsamplesinchain(Val{:burnin},10,20,5) == (50,50)
@test GeneralizedMetropolisHastings._numsamplesinchain(Val{:main},10,20,5) == (100,0)
@test GeneralizedMetropolisHastings._numsamplesinchain(Val{:all},10,20,5) == (151,51)

mhpolicy1 = policy(:mh,1)
mhpolicy2 = policy(:mh,1,store=:all)
mhpolicy3 = policy(:mh,1,store=:burnin,initialize=:default)

@test GeneralizedMetropolisHastings._storeduring(:burnin,mhpolicy1) == false
@test GeneralizedMetropolisHastings._storeduring(:burnin,mhpolicy2) == true
@test GeneralizedMetropolisHastings._storeduring(:burnin,mhpolicy3) == true
@test GeneralizedMetropolisHastings._storeduring(:main,mhpolicy1) == true
@test GeneralizedMetropolisHastings._storeduring(:main,mhpolicy2) == true
@test GeneralizedMetropolisHastings._storeduring(:main,mhpolicy3) == false

#create samplerstates
mhsamplerstate1 = samplerstate(rsampler1,smhnprops1,Float64,Float64)
mhsamplerstate2 = samplerstate(rsampler1,gmhnprops1,Float64,Float64)

#create all required objects
mhrunner1 = GMHRunnersTest.TestMHRunner(rnburnin1,rniter1,mhpolicy1)
mhrunner2 = GMHRunnersTest.TestMHRunner(rnburnin1,rniter1,mhpolicy2)
mhrunner3 = GMHRunnersTest.TestMHRunner(rnburnin1,rniter1,mhpolicy3)
mhindicator1,mhtunerstate1,mhchain1 = GeneralizedMetropolisHastings.createcommon(mhrunner1,rtuner1,rnparas1,smhnprops1,smhnprops1)
mhindicator2,mhtunerstate2,mhchain2 = GeneralizedMetropolisHastings.createcommon(mhrunner2,rtuner2,rnparas1,gmhnprops1,gmhnprops1)

@test numproposals(mhindicator1) == smhnprops1 && numsamples(mhindicator1) == smhnprops1
@test numproposals(mhindicator2) == gmhnprops1 && numsamples(mhindicator2) == gmhnprops1
@test indicatorsamples(mhindicator1) == [0,0]
@test indicatorsamples(mhindicator2) == [0,0,0]
@test numtunesteps(mhtunerstate1) == 2
@test numtunesteps(mhtunerstate2) == 1
@test numparas(mhchain1) == rnparas1 && numsamples(mhchain1) == rniter1
@test numparas(mhchain2) == rnparas1 && numsamples(mhchain2) == (rnburnin1 + rniter1)*gmhnprops1 + 1

#test the initialize function
#prepare some prior values
srand(rinitseed1) ; rpriorvals1 = initvalues(trait(:initialize,:prior),rparas1,mhpolicy1.sampletype)
srand(rinitseed1+1) ; rpriorvals2 = initvalues(trait(:initialize,:prior),rparas1,mhpolicy1.sampletype)
#copy (for testing) the first values of the chain
mcv1 = [mhchain1.values[1]]
mcll1 = mhchain1.loglikelihood[1]

#initialilze but without storing the first value in the chain
srand(rinitseed1) ; initialize!(mhrunner1,rmodel1,mhsamplerstate1,mhchain1,false)
mhll1 = [loglikelihood(rmodel1,evaluate!(rmodel1,rpriorvals1))]
GMHRunnersTest.testproposals(mhsamplerstate1,rpriorvals1,mhll1) #the value is firt generated in the proposals
GMHRunnersTest.testfrom(mhsamplerstate1,rpriorvals1,mhll1) #then copied to the from field
@test_approx_eq mean(mhsamplerstate1.density.distribution) rpriorvals1 #and the normal distribution of the mh sampler should also have been conditioned
GMHRunnersTest.teststore(mhchain1,mcv1,mcll1,1,0,0) #finally, the value has not been stored inside the chain

srand(rinitseed1+1) ; initialize!(mhrunner1,rmodel1,mhsamplerstate1,mhchain1,true)
mhll2 = [loglikelihood(rmodel1,evaluate!(rmodel1,rpriorvals2))]
GMHRunnersTest.testproposals(mhsamplerstate1,rpriorvals2,mhll2) #same tests as above
GMHRunnersTest.testfrom(mhsamplerstate1,rpriorvals2,mhll2)
@test_approx_eq mean(mhsamplerstate1.density.distribution) rpriorvals2 #and the normal distribution of the mh sampler should also have been conditioned
GMHRunnersTest.teststore(mhchain1,mhsamplerstate1.from.values,mhsamplerstate1.from.loglikelihood[1],1,0,1) #this time, the value should have been copied to the chain
