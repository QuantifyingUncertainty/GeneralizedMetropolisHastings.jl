###Test ODEModel using a spring/mass dynamic system
t1 = collect(0:0.1:10)
i1 = [-1.0,1.0]
p1 = [100.0,10.0]
n1 = [1e-1 1e-1]

m1 = spring_mass_model(t1,i1,p1,n1)
r1 = evaluate(m1,p1)

###Test the evaluate function for ODEModel
@test_approx_eq_eps evaluate(m1,p1) m1.measurements 1.0
@test_approx_eq_eps evaluate(m1,p1) spring_mass_data(t1,i1,[110.0,11.0],n1) 1.0
@test_approx_eq_eps evaluate(m1,[110.0,11.0]) spring_mass_data(t1,i1,[110.0,11.0],n1) 1.0

###Test the loglikelihood function for ODEModel
@test_approx_eq_eps loglikelihood(m1,r1)/length(r1) 0.0 0.2
@test_approx_eq_eps loglikelihood(m1,evaluate(m1,[112.0,11.0]))/length(r1) 0.0 0.4

###Test the logprior! function (for MCModels, but testing it with ODEModel)
s1 = BaseSample(p1)
@test ~isfinite(s1.logprior)
logprior!(m1,s1)
@test isfinite(s1.logprior) && s1.logprior == logprior!(m1,BaseSample(initvalues(ValuesFromPrior(),m1.parameters)))
@test isfinite(logprior!(m1,BaseSample([80.0,10.0]))) && ~isfinite(logprior!(m1,BaseSample([79.0,10.0])))

###Test the loglikelihood! function (for MCModels, but testing it with ODEModel)
@test ~isfinite(s1.loglikelihood)
r2 = loglikelihood!(m1,s1)
@test isfinite(s1.loglikelihood) && r2 == r1 && s1.loglikelihood == loglikelihood(m1,r1)

###Test the update_geometry function (for MCModels and BaseSample)
s2 = BaseSample(p1)
update_geometry!(m1,s2)
@test s1 == s2

show(s1)

####################

###Test the gradloglikelihood function for ODEModel
r11 = (evaluate(m1,[p1[1]+m1.gradientepsilon,p1[2]]) - r1)/m1.gradientepsilon
r12 = (evaluate(m1,[p1[1],p1[2]+m1.gradientepsilon]) - r1)/m1.gradientepsilon
@test_approx_eq_eps gradloglikelihood(m1,r1,cat(3,r11,r12))/length(r1) zeros(2) 0.2

###Test the gradlogprior! function (for MCModels, but testing with ODEModel)
g1 = GradientSample(p1)
@test all(~isfinite(g1.gradlogprior))
gradlogprior!(m1,g1)
@test g1.gradlogprior == zeros(size(g1.gradlogprior)) #uniform priors do not introduce a gradient
g1 = GradientSample([120.0,10.0])
@test all(~isfinite(g1.gradlogprior))
gradlogprior!(m1,g1)
@test ~isfinite(g1.gradlogprior[1]) && isfinite(g1.gradlogprior[2]) #parameter value is on the boundary of the prior, should become infinite

###Test the gradloglikelihood! function (for MCModels, but testing it with ODEModel)
g1 = GradientSample(p1)
logprior!(m1,g1)
gradlogprior!(m1,g1)
r1 = loglikelihood!(m1,g1)
f1 = gradient!(m1,g1,r1)
@test_approx_eq_eps f1 cat(3,r11,r12) 1e-6
g2 = GradientSample(p1)
update_geometry!(m1,g2)
@test g1 == g2

show(g1)

######################

###Test the tensorvalue function for the ODEModel
tv1 = tensorvalue(m1,cat(3,r11,r12),1,1)
@test_approx_eq_eps tv1/length(r1) 0.0 0.4

###Test the tensorlogprior! function (for MCModels, but testing with ODEModel)
ts1 = TensorSample(p1)
@test tensorlogprior!(m1,ts1) == zeros(2,2)

###Test the tensor! function (for MCModels, but testing with ODEModel)
ts1 = TensorSample(p1)
logprior!(m1,ts1)
gradlogprior!(m1,ts1)
tensorlogprior!(m1,ts1)
r1 = loglikelihood!(m1,ts1)
f1 = gradient!(m1,ts1,r1)
tensor!(m1,ts1,r1,f1)
@test_approx_eq_eps ts1.tensorloglikelihood[1,1] tv1 1e-6
ts2 = TensorSample(p1)
update_geometry!(m1,ts2)
@test ts1 == ts2

show(ts1)

###Test the approximte tensor calculation for the ODEModel
tg1 = GeneralizedMetropolisHastings.tangentvector(m1,r1,f1)
@test length(tg1) == 2

###Test the tensor! function for approximate tensor calculation (for MCModels, but testing with ODEModel)
srand(715)
a1 = ApproximateTensorSample(p1,100)
logprior!(m1,a1)
gradlogprior!(m1,a1)
tensorlogprior!(m1,a1)
r1 = loglikelihood!(m1,a1)
f1 = gradient!(m1,a1,r1)
tensor!(m1,a1,r1,f1)
@test size(a1.tangentvectors) == (2,100)
srand(715)
a2 = ApproximateTensorSample(p1,100)
update_geometry!(m1,a2)
@test_approx_eq_eps (ts1.tensorloglikelihood-a1.tensorloglikelihood)./ts1.tensorloglikelihood zeros(2,2) 0.6 #TODO: approximation is only half of the exact value. which one is right?
@test_approx_eq_eps (ts1.tensorloglikelihood-a2.tensorloglikelihood)./ts1.tensorloglikelihood zeros(2,2) 0.6

show(a1)















