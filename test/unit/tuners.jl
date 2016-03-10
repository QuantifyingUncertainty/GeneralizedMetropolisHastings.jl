t1 = tuner(:monitor,10)
t2 = tuner(:scale,4,0.7,:logistic)
t3 = tuner(:scale,4,0.5,:erf)

@test period(t1) == 10
@test period(t2) == 4
@test period(t3) == 4

@test verbose(t1)
@test verbose(t2)
@test verbose(t3)

s1 = tunerstate(t1)
s2 = tunerstate(t2)
s3 = tunerstate(t3)

@test isnan(rate(s1))
@test isnan(rate(s2))
@test isnan(rate(s3))

s1.proposed = 10 ; s1.accepted = 8 ; s1.totalproposed = 20
s2.proposed = 40 ; s2.accepted = 25 ; s2.totalproposed = 80
s3.proposed = 20 ; s3.accepted = 11 ; s3.totalproposed = 40

@test proposed(s1) == 10 && accepted(s1) == 8 && rate(s1) == 8/10 && totalproposed(s1) == 20
@test proposed(s2) == 40 && accepted(s2) == 25 && rate(s2) == 25/40 && totalproposed(s2) == 80
@test proposed(s3) == 20 && accepted(s3) == 11 && rate(s3) == 11/20 && totalproposed(s3) == 40

@test tune(t1,s1) == ()
@test tune(t2,s2) < 1.0
@test tune(t3,s3) > 1.0
@test_approx_eq tune(t2,s2) GeneralizedMetropolisHastings.tunelogistic(accepted(s2)/proposed(s2)-t2.targetrate,7.0)
@test_approx_eq tune(t3,s3) GeneralizedMetropolisHastings.tuneerf(accepted(s3)/proposed(s3)-t3.targetrate,3.0)

accepted!(s1,IndicatorStationary{Float64}(zeros(6),[1,1,1,2,3,4]))
accepted!(s2,IndicatorStationary{Float64}(zeros(5),[1,2,3,4,5]))
accepted!(s3,IndicatorStationary{Float64}(zeros(11),ones(Int,11)))

@test proposed(s1) == 15 && accepted(s1) == 11 && totalproposed(s1) == 25
@test proposed(s2) == 44 && accepted(s2) == 29 && totalproposed(s2) == 84
@test proposed(s3) == 30 && accepted(s3) == 11 && totalproposed(s3) == 50

resetburnin!(s1)
resetburnin!(s2)
resetburnin!(s3)

@test proposed(s1) == 0 && accepted(s1) == 0 && totalproposed(s1) == 25
@test proposed(s2) == 0 && accepted(s2) == 0 && totalproposed(s2) == 84
@test proposed(s3) == 0 && accepted(s3) == 0 && totalproposed(s3) == 50

###test the show functions
println()
println("====================")
println("Test show() function")
println("====================")
show(t1)
show(t2)
show(t3)
show(s1)
show(s2)
show(s3)
println("====================")
println("End  show() function")
println("====================")
println()
