using GeneralizedMetropolisHastings

###tests of BaseSample
l1 = 2
v1 = [1.0,2.0]
n1 = fill(NaN,l1)
m1 = fill(NaN,l1,l1)

s1 = BaseSample(v1)
s2 = BaseSample(l1)

@test isequal(s1.values,v1)
@test isequal(s2.values,zeros(2))
@test s1 != s2
@test s1 == BaseSample(v1)
@test s2 == BaseSample(l1)

###tests of GradientSample
s3 = GradientSample(v1)
s4 = GradientSample(l1)

@test isequal(s3.values,v1) && isequal(s3.gradloglikelihood,n1) && isequal(s3.gradlogprior,n1)
@test isequal(s4.values,zeros(2)) && isequal(s4.gradloglikelihood,n1) && isequal(s4.gradlogprior,n1)
@test s3 != s4
@test s3 == GradientSample(v1)
@test s4 == GradientSample(l1)

###tests of TensorSample
s5 = TensorSample(v1)
s6 = TensorSample(l1)

@test isequal(s5.values,v1) && isequal(s5.gradloglikelihood,n1) && isequal(s5.gradlogprior,n1) && isequal(s5.tensorloglikelihood,m1) && isequal(s5.tensorlogprior,m1)
@test isequal(s6.values,zeros(2)) && isequal(s6.gradloglikelihood,n1) && isequal(s6.gradlogprior,n1) && isequal(s6.tensorloglikelihood,m1) && isequal(s6.tensorlogprior,m1)
@test s5 != s6
@test s5 == TensorSample(v1)
@test s6 == TensorSample(l1)

###tests of vector equality
@test fill(BaseSample(2),2) == fill(BaseSample(2),2)
@test fill(BaseSample(2),2) != fill(BaseSample(3),2)

println()
println("====================")
println("Test show() function")
show(s1)
show(s3)
show(s5)
show([s1,s2])
println("End  show() function")
println("====================")
println()

nothing
