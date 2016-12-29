#specific multi-process tests of random generators, jobsegments and runners
@test nprocs() > 1 #only run this test with more than 1 process

#get the first worker process
where1 = workers()[1]
println("===============================================")
println("Multi-process tests with ",nprocs()," processes")
println("===============================================")
println("Local process: 1")
println("Remote process: ",where1)
println("========================================")
println("Multi-process tests of random generators")
println("========================================")
#but when the seed has been set in both, they should be the same
rng1loc = srand(345)
rng1rem = remotecall_fetch(srand,where1,345)
@test rng1loc.seed == rng1rem.seed
v1loc = rand()
v1rem = remotecall_fetch(rand,where1)
@test v1loc == v1rem

#each random generator should now be in the same state
v2loc = rand()
v2rem = remotecall_fetch(rand,where1)
@test v2loc == v2rem

#now bring them into different states
rng3loc = srand()
rng3rem = remotecall_fetch(srand,where1)
@test rng3loc.seed != rng3rem.seed
v3loc = rand()
v3rem = remotecall_fetch(rand,where1)
@test v3loc != v3rem

################################################################################

println("==============================================")
println("Multi-process tests of sampler state functions")
println("==============================================")
#create a sampler and sampler states
s1loc = sampler(:mh,:normal,0.1,2)
#create an indicatorstate
ind1loc = samplerstate(s1loc,1,Float64,Float64,false)
#create a local and remote auxiliary state
aux1loc = samplerstate(s1loc,3,Float64,Float64,true)
aux1rem = samplerstate(s1loc,3,Float64,Float64,true)

#set the local and remote random generator in the same state
rng4loc = srand(789)
rng4rem = remotecall_fetch(srand,where1,789)
@test rng4loc.seed == rng4rem.seed

#test that the propose! function gives the same results on local and remote
aux2loc = propose!(aux1loc)
aux2rem = remotecall_fetch(propose!,where1,aux1rem)
@test aux2loc.from.values == aux2rem.from.values
@test aux2loc.proposals.values == aux2rem.proposals.values
aux3loc = propose!(aux2loc)
aux3rem = remotecall_fetch(propose!,where1,aux2rem)
@test aux3loc.from.values == aux3rem.from.values
@test aux3loc.proposals.values == aux3rem.proposals.values

#set the from field of auxiliary states via indicator state
propose!(ind1loc)
aux4loc = prepareauxiliary!(ind1loc,aux3loc)
aux4rem = remotecall_fetch(prepareauxiliary!,where1,ind1loc,aux3rem)
@test aux4loc.from.values == ind1loc.proposals.values
@test aux4rem.from.values == ind1loc.proposals.values
@test aux4loc.proposals.values == aux4rem.proposals.values

#repeat propose!
rng5loc = srand(78932)
rng5rem = remotecall_fetch(srand,where1,78932)
aux5loc = propose!(aux4loc)
aux5rem = remotecall_fetch(propose!,where1,aux4rem)
@test aux5loc.from.values == aux5rem.from.values
@test aux5loc.proposals.values == aux5rem.proposals.values

#now set the local and remote random generators to different seeds
rng6loc = srand()
rng6rem = remotecall_fetch(srand,where1)
aux6loc = propose!(aux5loc)
aux6rem = remotecall_fetch(propose!,where1,aux5rem)
@test aux6loc.from.values == aux6rem.from.values
@test aux6loc.proposals.values != aux6rem.proposals.values

################################################################################
println("================================================")
println("Multi-process tests of remote segement functions")
println("================================================")
########Create all required objects
nprops1 = 5
niter1 = 6
nburnin1 = 4
ntunerperiod = 2
time1 = linspace(0.0,10.0,100)
params1 = parameters([:a,:b],[Distributions.Normal(3.0,0.2),Distributions.Normal(4.0,0.2)],[3.0,4.0])
nparams1 = length(params1)
data1 = data(:array,time1,cat(2,sin(3time1),cos(4time1)))
noise1 = noise(:gaussian,[0.01,0.01])
model1 = model(:target,params1,data1,noise1,(p,t)->cat(2,sin(p[1]*t),cos(p[2]*t)),time1;name="MultiProcessTestModel")
policy1 = policy(:mh,nprops1;jobsegments=:one,store=:all)
runner1 = runner(policy1,niter1,nprops1;numburnin =nburnin1)
@test GeneralizedMetropolisHastings._numjobsegments(policy1,nprops1) == 1
@test collect(GeneralizedMetropolisHastings._processnumbers(policy1,1)) == [where1]
tuner1 = tuner(:scale,ntunerperiod,0.5,:erf)
ind1loc = samplerstate(s1loc,1,Float64,Float64,false)
ind1rem = samplerstate(s1loc,1,Float64,Float64,false)
seg1loc = segment(policy1,model1,s1loc,nprops1)
seg1rem = remotesegments(policy1,model1,s1loc,nprops1)
indicator1loc,ts1loc,chain1loc = GeneralizedMetropolisHastings.createcommon(runner1,tuner1,nparams1,nprops1,nprops1)
indicator1rem,ts1rem,chain1rem = GeneralizedMetropolisHastings.createcommon(runner1,tuner1,nparams1,nprops1,nprops1)

#initialize indicator state and generate auxiliary variable
srand(383928) ; initialize!(runner1,model1,ind1loc,chain1loc,true)
srand(383928) ; initialize!(runner1,model1,ind1rem,chain1rem,true)

function testoneiteration(indloc,indrem,segloc,segrem,indicatorloc,indicatorrem,r,m,c,s1,s2,w,tl,tr)
    srand(s1) ; auxiliary!(r,m,indloc) ; indaccloc = acceptance!(indloc)
    srand(s1) ; auxiliary!(r,m,indrem) ; indaccrem = acceptance!(indrem)
    @test from(indloc) == from(indrem)
    @test proposals(indloc) == proposals(indrem)
    @test indaccloc == indaccrem
    srand(s2) ; accloc = iterate!(segloc,indloc)
    remotecall_wait(srand,w,s2) ; accrem = iterate!(segrem,indrem)
    sleep(1.0)
    @test isequal(accloc,accrem)
    transitionprobability!(indicatorloc,indaccloc,accloc) ; sampleindicator!(indicatorloc)
    transitionprobability!(indicatorrem,indaccrem,accrem) ; sampleindicator!(indicatorrem) ; store!(r,indrem,segrem,indicatorrem,c)
    @test indicatorloc.stationary == indicatorrem.stationary
    @test indicatorloc.samples == indicatorrem.samples
    ###Print
    println("Indicator state from: ",from(indrem).values)
    println("Indicator state proposals: ",proposals(indrem).values)
    println("Indicator state acceptance: ",indaccrem)
    println("Segment state acceptance: ",accrem)
    println("Stationary indicator: ",indicatorrem.stationary)
    println("Indicator samples: ",indicatorrem.samples)
    println("Remote segment acceptances: ",segrem.acceptances)
    println("Remote segment collected samples: ")
    dump(segrem.collectedsamples)
    println("Remote segment prop2collected: ")
    show(segrem.prop2collected)
    println()
    ###Prepare for next
    accepted!(tl,indicatorloc) ; indlocnew = preparenext!(r,indloc,segrem,indicatorloc)
    accepted!(tr,indicatorrem) ; indremnew = preparenext!(r,indrem,segrem,indicatorrem)
    @test from(indlocnew) == from(indremnew)
    @test proposals(indlocnew) == proposals(indremnew)
    println("Indicator state from new: ",from(indremnew).values)
    println("Indicator state proposals new: ",proposals(indremnew).values)
    indlocnew,indremnew
end

for i=1:10
    indloc,indrem = testoneiteration(ind1loc,ind1rem,seg1loc,seg1rem,indicator1loc,indicator1rem,runner1,model1,chain1rem,750-i,190+i,where1,ts1loc,ts1rem)
end

println("==========================")
println("Multi-process tests passed")
println("==========================")
nothing
