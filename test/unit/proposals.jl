import GeneralizedMetropolisHastings.update_density!
import GeneralizedMetropolisHastings.propose!
import GeneralizedMetropolisHastings.propose

m1 = zeros(2)
m2 = ones(2)
Σ1 = eye(2)
Σ2 = [1.0 0.0;0.0 4.0]
Σ3 = 4.0*eye(2)

###test the constructors
n1 = NormalDensity(Σ1)

@test length(n1.normal) == 2 && mean(n1.normal) == m1 && cov(n1.normal) == Σ1

n2 = NormalDensity(m2,Σ2)
n3 = NormalDensity(m1,Σ3)

@test length(n2.normal) == 2 && mean(n2.normal) == m2 && cov(n2.normal) == Σ2

###test the equality operators
@test NormalDensity(eye(3)) == NormalDensity(zeros(3),eye(3))
@test NormalDensity(eye(3)) != NormalDensity(ones(3),eye(3))
@test NormalDensity(eye(3)) != NormalDensity(eye(2))

###test the update functions
@test n1 != n2 != n3
update_density!(n1,m2)
update_density!(n2,Σ1)
update_density!(n3,m2,Σ1)
@test n1 == n2 == n3

###test the propose that returns and the in-place propose functions
srand(0); s1 = BaseSample(propose(n1))
srand(0); s2 = BaseSample(2); propose!(n1,s2); s2

@test s1.values == s2.values
@test size(propose(n1,4)) == (2,4) #generate multiple values at once

###test the vectorized propose
s3 = [BaseSample(2) for i=1:2]
@test s3[1].values == s3[2].values == zeros(2)
propose!(n1,s3)
@test s3[1].values != s3[2].values
n4 = [NormalDensity(2) for i=1:2]
println(typeof(n4))
propose!(n4,s3)
@test s3[1].values != s3[2].values

###test the show functions
println()
println("====================")
println("Test show() function")
show(n1)
show(n4)
println("End  show() function")
println("====================")
println()

###testing the sampling behaviour
imax = 1000
m3 = [4.0,2.0,0.0]
n5 = NormalDensity(m3,eye(3))
s4 = [BaseSample(3) for i=1:imax]
propose!(n5,s4)
@test_approx_eq_eps mapreduce((x)->(x.values),+,s4)/imax m3 0.1

###test the calculcation of proposal probabilities
n10 = NormalDensity([0.0,0.0],eye(2))
n11 = NormalDensity([2.0,2.0],eye(2))
s10 = BaseSample([0.0,0.0])
s11 = BaseSample([2.0,2.0])
@test logprobability(n10,s11) == logprobability(n11,s10)
@test logprobability(n10,s10) != logprobability(n10,s11)
@test logprobability([n10,n11],s11) == logprobability([n11,n10],s10)

nothing
