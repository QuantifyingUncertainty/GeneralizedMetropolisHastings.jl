#Define a very simple TargetModel for all tests below
t1 = linspace(0.0,10.0,100)
p1 = parameters([:a],[1.0],[5.0],[3.0])
d1 = data(:array,t1,sin(3t1))
n1 = noise(:gaussian,[0.01])
m1 = model(:target,p1,d1,n1,(p,t)->sin(p[1]t),t1;name="JobSegmentsTestModel")

#Define an MHSampler with normal proposal distribution
s1 = sampler(:mh,:normal,0.1,eye(1))

#Test a single segment with 2 proposals
println("==================")
println("Testing GMHSegment")
println("==================")
nprops1 = 2
nparas1 = length(p1)
policy1 = policy(:mh,nprops1)
indicatorstate1 = samplerstate(s1,1,policy1.sampletype,policy1.calculationtype)
auxiliarystate1 = samplerstate(s1,nprops1,policy1.sampletype,policy1.calculationtype)
copy!(indicatorstate1.proposals.values,[3.0])

segment1 = segment(policy1,m1,s1,nprops1)

@test numproposals(segment1) == 2
@test segment1.samplerstate.from.values == [0.0]
srand(234) ; prepareauxiliary!(indicatorstate1,auxiliarystate1) ; propose!(auxiliarystate1) ; geometry!(m1,proposals(auxiliarystate1)) ; a1 = acceptance!(auxiliarystate1)
srand(234) ; a2 = iterate!(segment1,indicatorstate1)
@test segment1.samplerstate.from.values == [3.0]
@test a1 == a2

@test getsamples(segment1,[2,1,2]) == copy(proposals(auxiliarystate1),[2,1,2])
@test getsamples(segment1,2) == copy(proposals(auxiliarystate1),2)
@test getsamples(segment1,[]) == samples(:base,1,0,Float64,Float64)
@test getsamplerstatevars(segment1)["density"] == density(segment1.samplerstate)

@test_approx_eq cov(segment1.samplerstate.density.distribution) cov(auxiliarystate1.density.distribution)
@test_approx_eq cov(segment1.samplerstate.density.distribution) 0.01*eye(1)
tune!(segment1,2.0)
tune!(auxiliarystate1,2.0)
@test_approx_eq cov(segment1.samplerstate.density.distribution) cov(auxiliarystate1.density.distribution)
@test_approx_eq cov(segment1.samplerstate.density.distribution) 0.04*eye(1)

println()
println("====================")
println("Test show() function")
show(segment1)
println("End  show() function")
println("====================")
println()

println("======================")
println("Testing RemoteSegments")
println("======================")

#test of number of job segments
@test GeneralizedMetropolisHastings._numjobsegments(Val{:procs}) == nprocs()
@test GeneralizedMetropolisHastings._numjobsegments(Val{:workers}) == nworkers()
@test GeneralizedMetropolisHastings._numjobsegments(Val{:one}) == 1
@test GeneralizedMetropolisHastings._numjobsegments(Val{:two}) == 2

#test the process numbers on which job segments will run
@test GeneralizedMetropolisHastings._processnumbers(Val{:procs}) == procs()
@test GeneralizedMetropolisHastings._processnumbers(Val{:workers}) == workers()
@test GeneralizedMetropolisHastings._processnumbers(Val{:one}) == workers()
@test GeneralizedMetropolisHastings._processnumbers(Val{:two}) == workers()

#test how many proposals per segment
nproposals = [2,3,4,5,6,7]
nsegments = [2,3,3,3,3,3]
npropsperseg = [1,1,2,2,2,3]
for i=1:length(nproposals)
    @test GeneralizedMetropolisHastings._numproposalspersegment(nproposals[i],nsegments[i]) == npropsperseg[i]
end

rprops1 = 5
rpolicy1 = policy(:mh,rprops1,jobsegments=:two)
rsegments1 = remotesegments(rpolicy1,m1,s1,rprops1)
@test numsegments(rsegments1) == 2
@test numproposalspersegment(rsegments1) == 3
@test numtotalproposals(rsegments1) == 6
@test GeneralizedMetropolisHastings._numjobsegments(rpolicy1,rprops1) == 2
@test collect(GeneralizedMetropolisHastings._processnumbers(rpolicy1,2)) == GMHRunnersTest.gettwoprocessnumbers(nprocs())
@test GeneralizedMetropolisHastings._insegmentindex(rsegments1,[1,4,5,6]) == Array{Int,1}[[1],[1,2,3]] #create the in-segment index for overall proposal index
@test GeneralizedMetropolisHastings._insegmentindex(rsegments1,[1,4,5,6,1,3,4,5,5]) == Array{Int,1}[[1,3],[1,2,3]] #repetitions are filtered out
@test rsegments1.prop2collected == Dict{Int,Tuple{Int,Int}}()
GeneralizedMetropolisHastings._insegmentindex2collected!(rsegments1,GeneralizedMetropolisHastings._insegmentindex(rsegments1,[1,4,5,6]))
@test rsegments1.prop2collected == Dict(1=>(1,1),4=>(1,2),5=>(2,2),6=>(3,2))

#prepare for following tests
indicatorstate1 = samplerstate(s1,1,policy1.sampletype,policy1.calculationtype)
initialize!(trait(:initialize,:default),m1,indicatorstate1.proposals)
geometry!(m1,indicatorstate1.proposals)
lsegments1 = [segment(policy1,m1,s1,numproposalspersegment(rsegments1)) for i=1:numsegments(rsegments1)]

#test the iterate function
srand(435) ; a1 = iterate!(rsegments1,indicatorstate1)
srand(435) ; a2 = cat(1,iterate!(lsegments1[1],indicatorstate1),iterate!(lsegments1[2],indicatorstate1))
@test nprocs() > 1 || a1 == a2 #only test this if all segments run on the same process

#test the prepare! function which copies a sample back into the indicator state
prepare!(rsegments1,indicatorstate1,2)
@test  nprocs() > 1 || indicatorstate1.from == copy(lsegments1[1].samplerstate.proposals,2)
prepare!(rsegments1,indicatorstate1,6)
@test  nprocs() > 1 || indicatorstate1.from == copy(lsegments1[2].samplerstate.proposals,3)

#test the retrievesamples! function which accesses the remote segments to retrieve samples from
samples1 = retrievesamples!(rsegments1,[1,4,5,6])
samples2 = map(getsamples,lsegments1,Array{Int,1}[[1],[1,2,3]])
@test  nprocs() > 1 || samples1 == samples2

#test the getsamples function for an index previously retrieved
@test  nprocs() > 1 || getsamples(rsegments1,4) == getsamples(lsegments1[2],1)

#test the getsamples function for an index not previously retrieved
@test  nprocs() > 1 || getsamples(rsegments1,3) == getsamples(lsegments1[1],3)

#test storing the smaples into the chain
c1 = chain(:standard,1,3,policy1.sampletype,policy1.calculationtype)
c2 = chain(:standard,1,3,policy1.sampletype,policy1.calculationtype)
store!(c1,indicatorstate1.proposals,1)
store!(c2,indicatorstate1.proposals,1)
store!(rsegments1,c1,1)
store!(c2,lsegments1[1].samplerstate.proposals,1)
@test  nprocs() > 1 || c1.values == c2.values

#test tuning the segments
tune!(rsegments1,0.5)
map((s)->tune!(s,0.5),lsegments1)
@test_approx_eq cov(lsegments1[1].samplerstate.density.distribution) cov(fetch(rsegments1.remote[1]).samplerstate.density.distribution)

if nprocs() == 1 #fix this later, currently fails for more than 1 process
    rseg1statevars = getsamplerstatevars(rsegments1)
    @test length(rseg1statevars) == 2
    @test isa(rseg1statevars[1]["density"],DistributionWrapper)
    @test isa(rseg1statevars[2]["density"],DistributionWrapper)
end

println()
println("====================")
println("Test show() function")
show(rsegments1)
nprocs()==1?showsamplerstatevars(rsegments1):nothing #fix this later
println("End  show() function")
println("====================")
println()

nothing
