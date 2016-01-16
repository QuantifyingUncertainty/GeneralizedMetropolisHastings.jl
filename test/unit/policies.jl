@test GeneralizedMetropolisHastings._num2proposefrom(1) == ProposeFromIndicator
@test GeneralizedMetropolisHastings._num2proposefrom(2) == ProposeFromAuxiliary
@test_throws AssertionError GeneralizedMetropolisHastings._num2proposefrom(0)

p1 = GeneralizedMetropolisHastings._policy(Val{:generic},InitializeFromDefault,ProposeFromAuxiliary,IndicatorStationary,Float64)
@test typeof(p1) <: GenericPolicy && p1.initialize == InitializeFromDefault && p1.propose == ProposeFromAuxiliary && p1.indicate == IndicatorStationary && p1.numbertype == Float64

p2 = policy(:generic,InitializeFromDefault,IndicatorStationary,10)
p3 = policy(:generic,InitializeFromPrior,IndicatorCyclical,1,Float32)

@test p1 == p2
@test typeof(p3) <: GenericPolicy && p3.initialize == InitializeFromPrior && p3.propose == ProposeFromIndicator && p3.indicate == IndicatorCyclical && p3.numbertype == Float32

println("====================")
println("Test show() function")
show(p1)
show(p3)
println("End  show() function")
println("====================")
println()

nothing
