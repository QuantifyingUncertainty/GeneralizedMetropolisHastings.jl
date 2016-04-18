t1 = tuner(:monitor,10)
t2 = tuner(:scale,4,0.7,:logistic)
t3 = tuner(:scale,4,0.5,:erf)

@test period(t1) == 10
@test period(t2) == 4
@test period(t3) == 4

@test verbose(t1)
@test verbose(t2)
@test verbose(t3)

s1 = tunerstate(t1,100,Float64)
s2 = tunerstate(t2,32,Float64)
s3 = tunerstate(t3,33,Float32)

@test all(isnan(rate(s1))) && numtunesteps(s1) == 10
@test all(isnan(rate(s2))) && numtunesteps(s2) == 8
@test all(isnan(rate(s3))) && numtunesteps(s3) == 9

s1.proposed[s1.index] = 10 ; s1.accepted[s1.index] = 8
s2.proposed[s2.index] = 40 ; s2.accepted[s2.index] = 25
s3.proposed[s3.index] = 20 ; s3.accepted[s3.index] = 11

@test rate(s1)[s1.index] == 8/10 && total(s1) == 10
@test rate(s2)[s2.index] == 25/40 && total(s2) == 40
@test rate(s3)[s3.index] == 11/20 && total(s3) == 20

@test tune(t1,s1) == ()
@test tune(t2,s2) < 1.0
@test tune(t3,s3) > 1.0
@test_approx_eq s2.scalefactor[s2.index] GeneralizedMetropolisHastings.tunelogistic(accepted(s2)[s2.index]/proposed(s2)[s2.index]-t2.targetrate,7.0)
@test_approx_eq s3.scalefactor[s3.index] GeneralizedMetropolisHastings.tuneerf(accepted(s3)[s3.index]/proposed(s3)[s3.index]-t3.targetrate,3.0)

nextindex!(s1)
nextindex!(s2)
nextindex!(s3)

@test index(s1) == 2 && index(s2) == 2 && index(s3) == 2

accepted!(s1,IndicatorStationary{Float64}(zeros(6),[1,1,1,2,3,4]))
accepted!(s2,IndicatorStationary{Float64}(zeros(5),[1,2,3,4,5]))
accepted!(s3,IndicatorStationary{Float64}(zeros(11),ones(Int,11)))

@test proposed(s1)[s1.index] == 5 && accepted(s1)[s1.index] == 3 && total(s1) == 15
@test proposed(s2)[s2.index] == 4 && accepted(s2)[s2.index] == 4 && total(s2) == 44
@test proposed(s3)[s3.index] == 10 && accepted(s3)[s3.index] == 0 && total(s3) == 30

s2.index = 9
println("A warning should be printed below this line")
accepted!(s2,IndicatorStationary{Float64}(zeros(5),[1,2,3,4,5]))

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
showstep(t1,s1)
showstep(t2,s2)
showstep(t3,s3)
println("====================")
println("End  show() function")
println("====================")
println()
