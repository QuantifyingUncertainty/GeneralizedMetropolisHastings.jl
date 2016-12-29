###Define functions to be used in testing below

###Test functions for TargetModel
function targetsin!(r::Vector,paras::Vector,t::AbstractVector)
    for i=1:length(t)
        r[i] = sin(paras[1]*t[i])
    end
    r
end

targetsin(paras::Vector,t::AbstractVector) = targetsin!(zeros(eltype(t),length(t)),paras,t)

###Test function for ODEModel
function odesin(t,y,ydot,paras)
    ydot[1] = paras[1]*cos(paras[1]*t)
end

#########################################################################################################

### Test model creation and evaluation functions

timepoints1 = linspace(0.0,10.0,100)
parameters1 = parameters([:a],[1.0],[5.0],[3.0])
measurements1 = data(:function,timepoints1,targetsin,[3.0],timepoints1)
measurements2 = data(:function,timepoints1,targetsin,[2.5],timepoints1)
measurements3 = data(:function,timepoints1,targetsin,[2.0],timepoints1)
noisemodel1 = noise(:gaussian,[0.01])
initial1 = [0.0]

model1 = model(:target,parameters1,measurements1,noisemodel1,targetsin,timepoints1;name="Test")
model2 = model(:target!,parameters1,measurements1,noisemodel1,targetsin!,timepoints1;name="Test!")
model3 = model(:ode,parameters1,measurements1,noisemodel1,odesin,initial1,1,[1];name="Test")

@test_approx_eq_eps evaluate!(model1,[3.0]) measurements(model1) 1e-4
@test_approx_eq_eps evaluate!(model2,[3.0]) measurements(model2) 1e-4
@test_approx_eq_eps evaluate!(model3,[3.0]) measurements(model3) 1e-4

#########################################################################################################

### Test the common initialize! function
s0 = samples(:base,1,1,Float64,Float64)
@test initialize!(trait(:initialize,:default),model1,s0).values == [3.0]
@test (srand(345) ; initialize!(trait(:initialize,:prior),model1,s0).values) == (srand(345) ; [rand(parameters1[1].prior)])

### Test the common interface functions
type TestModel <: AbstractModel
    parameters::Vector
end

@test numparas(TestModel(parameters1)) == 1
@test parameters(TestModel(parameters1)) == parameters1

@test_throws MethodError evaluate!(TestModel(parameters1),[1.0])
@test_throws MethodError dataindex(TestModel(parameters1))
@test_throws MethodError measurements(TestModel(parameters1))
@test_throws MethodError noisemodel(TestModel(parameters1))

########################################################################################################

### Test the geometry calculations
r1 = copy(evaluate!(model1,[3.0]))
r2 = copy(evaluate!(model2,[2.5]))
r3 = copy(evaluate!(model3,[2.0]))

@test_approx_eq loglikelihood(model1,r1) loglikelihood(noisemodel1,datavalues(measurements1),r1)
@test_approx_eq loglikelihood(model2,r2) loglikelihood(noisemodel1,datavalues(measurements1),r2)
@test_approx_eq loglikelihood(model3,r3) loglikelihood(noisemodel1,datavalues(measurements1),r3)

### Test the geometry function for zero-order samples
for m in [model1,model2,model3]
    s = samples(:base,1,3,Float64,Float64)
    copy!(s.values,[3.0,2.0,0.0])
    geometry!(m,s)
    ll = [loglikelihood(m,evaluate!(m,[3.0])) loglikelihood(m,evaluate!(m,[2.0])) -Inf]
    @test_approx_eq s.logprior logprior(parameters1,s.values,Float64)
    @test_approx_eq s.loglikelihood ll
end

###test the show functions
println()
println("====================")
println("Test show() function")
println("====================")
show(model1)
show(model2)
show(model3)
println("====================")
println("End  show() function")
println("====================")
println()
