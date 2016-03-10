###Test the trait types
for args in [(:initialize,InitializeFrom,[:default,:prior]),
            (:propose,ProposeFrom,[:indicator,:auxiliary]),
            (:indicator,IndicatorType,[:stationary,:cyclical]),
            (:samplerstates,SamplerStates,[:nprocs,:nworkers,:test])]
    for j in args[3]
        t = trait(args[1],j)
        @test typeof(t) == args[2]
        @test traitvalue(t) == j
        @test traittype(t) == Val{j}
    end
end

###Test GMHPolicy
@test traittype(GeneralizedMetropolisHastings._num2propose(1)) == Val{:indicator}
@test traittype(GeneralizedMetropolisHastings._num2propose(2)) == Val{:auxiliary}
@test_throws AssertionError GeneralizedMetropolisHastings._num2propose(0)

p1 = policy(:gmh,10)
p2 = policy(:gmh,1,initialize=:default,indicator=:cyclical,samplerstates=:test,sampletype=Int,calculationtype=Float32)

for args in [(p1,(:prior,:auxiliary,:stationary,:nprocs,Float64,Float64)),
             (p2,(:default,:indicator,:cyclical,:test,Int,Float32))]
    for (i,f) in enumerate([:initialize,:propose,:indicator,:samplerstates])
        @test traittype(getfield(args[1],f)) == Val{args[2][i]}
    end
    args[1].sampletype == args[2][5]
    args[1].calculationtype == args[2][6]
end

println("====================")
println("Test show() function")
show(p1)
show(p2)
println("End  show() function")
println("====================")
println()

nothing
