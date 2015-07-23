import GeneralizedMetropolisHastings.update_density!
import GeneralizedMetropolisHastings.propose!
import GeneralizedMetropolisHastings.propose

m1 = zeros(2)
m2 = ones(2)
Σ1 = eye(2)
Σ2 = [1.0 0.0;0.0 4.0]
Σ3 = 4.0*eye(2)
σ2 = [1.0,2.0]
s3 = 2.0

###test the constructors
n1 = NormalDensity(Σ1)
n2 = NormalDensity(σ2)
n3 = NormalDensity(2,s3)

@test length(n1.normal) == 2 && mean(n1.normal) == m1 && cov(n1.normal) == Σ1
@test length(n2.normal) == 2 && mean(n2.normal) == m1 && cov(n2.normal) == Σ2
@test length(n3.normal) == 2 && mean(n3.normal) == m1 && cov(n3.normal) == Σ3

n4 = NormalDensity(m2,Σ1)
n5 = NormalDensity(m2,σ2)
n6 = NormalDensity(m2,s3)

@test length(n4.normal) == 2 && mean(n4.normal) == m2 && cov(n4.normal) == Σ1
@test length(n5.normal) == 2 && mean(n5.normal) == m2 && cov(n5.normal) == Σ2
@test length(n6.normal) == 2 && mean(n6.normal) == m2 && cov(n6.normal) == Σ3

###test the equality operators
@test NormalDensity(eye(3)) == NormalDensity(zeros(3),eye(3))
@test NormalDensity(eye(3)) != NormalDensity(ones(3),eye(3))
@test NormalDensity(eye(3)) != NormalDensity(eye(2))

###test the update functions
update_density!(n1,m2,s3)
update_density!(n2,m2,σ2)
update_density!(n3,m2,Σ1)

@test length(n1.normal) == 2 && mean(n1.normal) == m2 && cov(n1.normal) == Σ3
@test length(n2.normal) == 2 && mean(n2.normal) == m2 && cov(n2.normal) == Σ2
@test length(n3.normal) == 2 && mean(n3.normal) == m2 && cov(n3.normal) == Σ1

update_density!(n4,s3)
update_density!(n5,σ2)
update_density!(n6,Σ1)

@test length(n4.normal) == 2 && mean(n4.normal) == m1 && cov(n4.normal) == Σ3
@test length(n5.normal) == 2 && mean(n5.normal) == m1 && cov(n5.normal) == Σ2
@test length(n6.normal) == 2 && mean(n6.normal) == m1 && cov(n6.normal) == Σ1

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
n7 = [NormalDensity(2) for i=1:2]
println(typeof(n7))
propose!(n7,s3)
@test s3[1].values != s3[2].values

###test the show functions
println()
println("====================")
println("Test show() function")
show(n1)
show(n7)
println("End  show() function")
println("====================")
println()
println("===================================================")
println("Test of info printout for 1-element constructors")
n7 = NormalDensity([0.0],[1.0])
n8 = NormalDensity([1.0])
println("End info printout. You should see 2 infos printed")
println("===================================================")

###testing the sampling behaviour
imax = 1000
m3 = [4.0,2.0,0.0]
n9 = NormalDensity(m3,eye(3))
s4 = [BaseSample(3) for i=1:imax]
propose!(n9,s4)
@test_approx_eq_eps mapreduce((x)->(x.values),+,s4)/imax m3 0.1

nothing
