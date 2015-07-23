import GeneralizedMetropolisHastings.MHHeap
import GeneralizedMetropolisHastings.SmMALAHeap
import GeneralizedMetropolisHastings.set_from!
import GeneralizedMetropolisHastings.propose!
import GeneralizedMetropolisHastings.update_proposal!

#############Test Metropolis-Hastings samplers####################

###test the MHNormal constructors
s1 = MHNormal(eye(2),2.0)
s2 = MHNormal(2,2.0)
@test s1.covariance == eye(2) && s1.initialscaling == 2.0
@test s1.covariance == s2.covariance && s1.initialscaling == s2.initialscaling

###test the MHHeap constructors and equality operators
h1 = MHHeap([BaseSample(2) for i=1:3],NormalDensity(2.0*eye(2)),zeros(2),2.0)
h2 = MHHeap(s1,3)
@test h1 ==  MHHeap([BaseSample(2) for i=1:3],NormalDensity(2.0*eye(2)),zeros(2),2.0)
@test h2 == MHHeap(s1,3)
@test h1 == h2

###test the from! functions
v1 = [1.0,2.0]
v2 = [4.0,0.0]
h1.samples[2].values = v1
b1 = BaseSample(v2)
set_from!(s1,h1,h1.samples[2])
@test h1.fromvalues == v1
set_from!(s1,h1,h1.samples[1])
@test h1.fromvalues == zeros(2)
set_from!(s1,h1,b1)
@test h1.fromvalues == v2

###test the propose! functions
@test h1.samples[1].values == zeros(2)
propose!(s1,h1,h1.samples[1])
@test h1.samples[1].values != zeros(2)

###test the update_proposal! function (should not have any effect)
h2 = deepcopy(h1)
update_proposal!(s1,h1,1)
@test h1 == h2

###test the show functions
println()
println("====================")
println("Test show() function")
show(h1)
println("End  show() function")
println("====================")
println()

imax = 10000
vals = Array(Float64,2,2*imax)
c = 0
for i=1:imax
  propose!(s1,h1,2)
  for j in [1,3]
    vals[:,c+=1] = h1.samples[j].values
  end
end
@test_approx_eq_eps mean(vals,2) v1 0.1
vals = Array(Float64,2,2*imax)
c = 0
for i=1:imax
  propose!(s1,h1,2,b1)
  for j in [1,3]
    vals[:,c+=1] = h1.samples[j].values
  end
end
@test_approx_eq_eps mean(vals,2) v2 0.1

###########################Test SmMALA samplers############################

###Test the sampler constructors
sm1 = SmMALANormal(2,0.1)
sm2 = TrSmMALANormal(2,0.1)
sm3 = TrSmMALARandomNormal(2,0.1)

@test typeof(sm1) <: GeneralizedMetropolisHastings.SmMALA && typeof(sm1) <: GeneralizedMetropolisHastings.SmMALANormalFamily
@test typeof(sm2) <: GeneralizedMetropolisHastings.TrSmMALA && typeof(sm2) <: GeneralizedMetropolisHastings.SmMALANormalFamily && typeof(sm2) <: GeneralizedMetropolisHastings.TrSmMALANormalFamily
@test typeof(sm3) <: GeneralizedMetropolisHastings.TrSmMALARandom && typeof(sm3) <: GeneralizedMetropolisHastings.SmMALANormalFamily && typeof(sm3) <: GeneralizedMetropolisHastings.TrSmMALANormalFamily
@test nparas(sm1) == nparas(sm2) == nparas(sm3) == 2

###Test the heap constructors
smh1 = SmMALAHeap(sm1,4)
smh2 = SmMALAHeap(sm2,4)
smh3 = SmMALAHeap(sm3,4)
@test smh1 == smh2 == smh3 #despite different samplers, all heaps should be the same

###Test the mean and covariance calculation functions
t1 = TensorSample(2)
t1.values = [0.0,0.1]
t1.loglikelihood = 0.3
t1.gradloglikelihood = [-0.5,0.5]
t1.tensorloglikelihood = [0.5 0.0;0.0 1.0]
m1 = [-0.005,0.1025]
c1 = [0.02 0.0;0.0 0.01]
m2 = [-0.047619,0.145455]
c2 = [0.0952381 0.0;0.0 0.0909091]
@test_approx_eq_eps GeneralizedMetropolisHastings.smmalamean(sm1,t1,0.1) m1 1.0e-6 #test of SmMALA mean
@test_approx_eq_eps GeneralizedMetropolisHastings.smmalacov(sm1,t1,0.1) c1 1.0e-6 #test of SmMALA cov
@test_approx_eq_eps GeneralizedMetropolisHastings.smmalamean(sm2,t1,0.1) m2 1.0e-6 #test of TrSmMALA mean
@test_approx_eq_eps GeneralizedMetropolisHastings.smmalacov(sm2,t1,0.1) c2 1.0e-6 #test of TrSmMALA cov
@test GeneralizedMetropolisHastings.smmalamean(sm2,t1,0.1) == GeneralizedMetropolisHastings.smmalamean(sm3,t1,0.1)
@test GeneralizedMetropolisHastings.smmalacov(sm2,t1,0.1) == GeneralizedMetropolisHastings.smmalacov(sm3,t1,0.1)

###Test the local update function
n1 = NormalDensity(2)
update_proposal!(sm1,n1,t1,0.1)
@test_approx_eq_eps mean(n1.normal) m1 1.0e-6 #test of SmMALA mean
@test_approx_eq_eps cov(n1.normal) c1 1.0e-6 #test of SmMALA cov
update_proposal!(sm2,n1,t1,0.1)
@test_approx_eq_eps mean(n1.normal) m2 1.0e-6 #test of TrSmMALA mean
@test_approx_eq_eps cov(n1.normal) c2 1.0e-6 #test of TrSmMALA cov

###Test the from! functions
set_from!(sm1,smh1,t1)
@test_approx_eq_eps mean(smh1.fromdensity.normal) m1 1.0e-6 #test of SmMALA mean
@test_approx_eq_eps cov(smh1.fromdensity.normal) c1 1.0e-6 #test of SmMALA cov
set_from!(sm2,smh2,t1)
@test_approx_eq_eps mean(smh2.fromdensity.normal) m2 1.0e-6 #test of SmMALA mean
@test_approx_eq_eps cov(smh2.fromdensity.normal) c2 1.0e-6 #test of SmMALA cov
set_from!(sm3,smh3,t1)
@test_approx_eq_eps mean(smh3.fromdensity.normal) m2 1.0e-6 #test of SmMALA mean
@test_approx_eq_eps cov(smh3.fromdensity.normal) c2 1.0e-6 #test of SmMALA cov

###Test the propose! functions
set_from!(sm1,smh1,t1)
@test_approx_eq_eps (srand(0); propose!(sm1,smh1,smh1.samples[1]); smh1.samples[1].values) (srand(0); rand(NormalDensity(m1,c1).normal)) 1.0e-6
set_from!(sm2,smh2,t1)
@test_approx_eq_eps (srand(0); propose!(sm2,smh2,smh2.samples[2]); smh2.samples[2].values) (srand(0); rand(NormalDensity(m2,c2).normal)) 1.0e-6
set_from!(sm3,smh3,t1)
@test (srand(0); propose!(sm2,smh2,smh2.samples[2]); smh2.samples[2].values == (srand(0); propose!(sm3,smh3,smh3.samples[3]); smh3.samples[3].values))

###test the show functions
println()
println("====================")
println("Test show() function")
show(smh1)
println("End  show() function")
println("====================")
println()

###Test the vectorized propose! functions
imax = 10000
vals = Array(Float64,2,3*imax)
c = 0
for i=1:imax
  propose!(sm1,smh1,2,t1)
  for j in [1,3,4]
    vals[:,c+=1] = smh1.samples[j].values
  end
end
@test_approx_eq_eps mean(vals,2) m1 0.1

imax = 10000
vals = Array(Float64,2,3*imax)
c = 0
for i=1:imax
  propose!(sm2,smh2,4,t1)
  for j in [1,2,3]
    vals[:,c+=1] = smh2.samples[j].values
  end
end
@test_approx_eq_eps mean(vals,2) m2 0.1
