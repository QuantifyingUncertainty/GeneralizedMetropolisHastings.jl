###Test the trait types
for args in [(:mhrunner,MHRunnerType,[:standard,:generalized]),
            (:initialize,InitializeFrom,[:default,:prior]),
            (:propose,ProposeFrom,[:indicator,:auxiliary]),
            (:indicator,IndicatorType,[:stationary,:cyclical]),
            (:jobsegments,JobSegments,[:procs,:workers,:none,:one,:two]),
            (:chain,ChainType,[:standard,:gradient]),
            (:store,StoreDuring,[:burnin,:main,:all])]
    for j in args[3]
        t = trait(args[1],j)
        @test typeof(t) == args[2]
        @test traitvalue(t) == j
        @test traittype(t) == Val{j}
    end
end

###Test MHRuntimePolicy
@test traitvalue(GeneralizedMetropolisHastings._num2propose(1)) == :indicator
@test traitvalue(GeneralizedMetropolisHastings._num2propose(2)) == :auxiliary
@test_throws AssertionError GeneralizedMetropolisHastings._num2propose(0)

@test traitvalue(GeneralizedMetropolisHastings._num2runner(1)) == :standard
@test traitvalue(GeneralizedMetropolisHastings._num2runner(2)) == :generalized
@test_throws AssertionError GeneralizedMetropolisHastings._num2runner(0)

@test traitvalue(GeneralizedMetropolisHastings._num2segments(1,:workers)) == :none
@test traitvalue(GeneralizedMetropolisHastings._num2segments(2,:workers)) == :workers
@test_throws AssertionError GeneralizedMetropolisHastings._num2segments(0,:none)

p1 = policy(:mh,10)
p2 = policy(:mh,1,initialize=:default,indicator=:cyclical,chain=:gradient,store=:all,sampletype=Int,calculationtype=Float32)

for args in [(p1,(:generalized,:prior,:auxiliary,:stationary,:workers,:standard,:main,Float64,Float64)),
             (p2,(:standard,:default,:indicator,:cyclical,:none,:gradient,:all,Int,Float32))]
    for (i,f) in enumerate([:runner,:initialize,:propose,:indicator,:jobsegments,:chain,:store])
        @test traittype(getfield(args[1],f)) == Val{args[2][i]}
    end
    @test args[1].sampletype == args[2][end-1]
    @test args[1].calculationtype == args[2][end]
end

println("====================")
println("Test show() function")
show(p1)
show(p2)
println("End  show() function")
println("====================")
println()

nothing
