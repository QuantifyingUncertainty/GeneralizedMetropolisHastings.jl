###The same model is used with the test cases below, a two-parameter spring-mass model
defaults1 = [110.0,8.0]
srand(897) ; model1 = springmassmodel(linspace(0.0,10.0,100),[-1.0,1.0],[100.0,10.0],[1e-2,9e-2],defaults1)
sampler1 = sampler(:mh,:normal,0.5,[1.0 0.0;0.0 0.1])
###################Test a runner with 1 proposal per iteration (normal Metropolis-Hastings)
println("==================")
println("Testing GMHSegment")
println("==================")
nprops1 = 2
nparas1 = 2
policy1 = policy(:mh,nprops1)
sampler1 = sampler(:mh,:normal,0.5,[1.0 0.0;0.0 0.1])
samplerstate1 = samplerstate(sampler1,nprops1,policy1.sampletype,policy1.calculationtype)
from1 = samples(:base,nparas1,1,policy1.sampletype,policy1.calculationtype) ; copy!(from1.values,defaults1)

segment1 = segment(policy1,model1,sampler1,nprops1)

@test numproposals(segment1) == 2

srand(234) ; setfrom!(samplerstate1,from1) ; propose!(samplerstate1) ; geometry!(model1,proposals(samplerstate1)) ; a1 = acceptanceratio!(samplerstate1)
srand(234) ; a2 = iterate!(segment1,from1)
@test a1 == a2

@test getsamples(segment1,[2,1,2]) == copy(proposals(samplerstate1),[2,1,2])
@test getsamples(segment1,2) == copy(proposals(samplerstate1),2)
@test getsamples(segment1,[]) == samples(:base,2,0,Float64,Float64)

tune!(segment1,sampler1,2.0)
tune!(sampler1,samplerstate1,2.0)

@test segment1.samplerstate.scalefactor == samplerstate1.scalefactor

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
@test GeneralizedMetropolisHastings._numjobsegments(Val{:test}) == 3

#test the process numbers on which job segments will run
@test GeneralizedMetropolisHastings._processnumbers(Val{:procs}) == procs()
@test GeneralizedMetropolisHastings._processnumbers(Val{:workers}) == workers()
@test GeneralizedMetropolisHastings._processnumbers(Val{:test}) == workers()

#define runtime variables
ntestrunners = 3
nproposals = [2,3,4]
nsegments = [2,3,3]
npropsperseg = [1,1,2]

testget = [2,3,4]
testindex = Vector{Vector{Int}}(3) ; testindex[1] = [2] ; testindex[2] = [1,3] ; testindex[3] = [1,3,4] ; testindex
testsegindex = Any[Vector{Vector{Int}}(2),Vector{Vector{Int}}(3),Vector{Vector{Int}}(3)]
testsegindex[1][1] = [] ; testsegindex[1][2] = [1]
testsegindex[2][1] = [1] ; testsegindex[2][2] = [] ; testsegindex[2][3] = [1]
testsegindex[3][1] = [1] ; testsegindex[3][2] = [1,2] ; testsegindex[3][3] = []
testcolindex = [Dict{Int,Tuple{Int,Int}}() for i=1:3]
testcolindex[1][2] = (1,2)
testcolindex[2][1] = (1,1) ; testcolindex[2][3] = (1,3)
testcolindex[3][1] = (1,1) ; testcolindex[3][3] = (1,2) ; testcolindex[3][4] = (2,2)
testsamples = Vector{Vector{Int}}(3) ; testsamples[1] = [2,2,2] ; testsamples[2] = [3,1,1,1,3] ; testsamples[3] = [1,3,4,4,1] ; testsamples

#create the runtime policies
policies = [policy(:mh,nproposals[i],jobsegments=:test) for i=1:ntestrunners]
chains = [chain(:standard,2,5,Float64,Float64) for i=1:ntestrunners]
segments = Vector{Any}(ntestrunners)

testprocs = repmat(workers(),GeneralizedMetropolisHastings._numjobsegments(Val{:test}))
indicatorstate1 = samplerstate(sampler1,1,Float64,Float64)
for i=1:ntestrunners
    #test creation aspects
    @test traitvalue(policies[i].runner) == :generalized
    @test GeneralizedMetropolisHastings._numjobsegments(policies[i],nproposals[i]) == nsegments[i]
    @test GeneralizedMetropolisHastings._numproposalspersegment(nproposals[i],nsegments[i]) == npropsperseg[i]
    procvals = collect(GeneralizedMetropolisHastings._processnumbers(policies[i],nsegments[i]))
    @test length(procvals) == nsegments[i] && procvals == testprocs[1:nsegments[i]]
    segments[i] = remotesegments(policies[i],model1,sampler1,nproposals[i])
    @test segments[i].numsegments == nsegments[i]
    for j=1:segments[i].numsegments
        @test isa(fetch(segments[i].remote[j]),GMHSegment)
        r = @spawnat segments[i].remote[j].where numproposals(fetch(segments[i].remote[j]))
        @test fetch(r) == npropsperseg[i]
    end
    #test iteration aspects
    a = iterate!(segments[i],proposals(indicatorstate1))
    @test length(a) == nsegments[i]
    for j=1:segments[i].numsegments
        @test length(a[i]) == npropsperseg[i]
    end
    sam1 = getsamples(segments[i],testget[i])
    #test getting the samples from the segments
    @test GeneralizedMetropolisHastings._prop2segment(segments[i],testindex[i]) == testsegindex[i]
    @test GeneralizedMetropolisHastings._segment2collected(segments[i],testsegindex[i]) == testcolindex[i]
    segments[i].prop2collected = Dict{Int,Tuple{Int,Int}}()
    b = retrievesamples!(segments[i],testindex[i])
    @test segments[i].prop2collected == testcolindex[i]
    for j=1:segments[i].numsegments
        @test numparas(b[j]) == 2 && numsamples(b[j]) == length(testsegindex[i][j])
    end
    p,s = testcolindex[i][testget[i]]
    @test sam1 == getsamples(segments[i],testget[i]) == copy(segments[i].collectedsamples[s],p)
    #test storing samples in the chain
    for j=1:length(testsamples[i])
        store!(segments[i],chains[i],testsamples[i][j])
    end
    for j=1:length(testsamples[i])
        p,s = segments[i].prop2collected[testsamples[i][j]]
        @test chains[i].values[:,j] == segments[i].collectedsamples[s].values[:,p]
    end
    tune!(segments[i],sampler1,1.0)
    for j=1:segments[i].numsegments
        @test fetch(segments[i].remote[j]).samplerstate.scalefactor == 0.5
    end
    tune!(segments[i],sampler1,1.5)
    for j=1:segments[i].numsegments
        @test fetch(segments[i].remote[j]).samplerstate.scalefactor == 0.75
    end
    tune!(segments[i],sampler1,0.5)
    for j=1:segments[i].numsegments
        @test fetch(segments[i].remote[j]).samplerstate.scalefactor == 0.375
    end
end

println()
println("====================")
println("Test show() function")
show(segments)
println("End  show() function")
println("====================")
println()

nothing
