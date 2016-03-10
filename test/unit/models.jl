###Test ODEModel using a spring/mass dynamic system
t1 = 0:0.1:10.0
i1 = [-1.0,1.0]
l1 = [80.0,8.0]
p1 = [100.0,10.0]
u1 = [120.0,12.0]
v1 = [1e-2,9e-2]

m1 = springmassmodel(t1,i1,p1,v1,l1,u1)
r1 = evaluate!(m1,p1)

@test numparas(m1) == 2

###Test the evaluate function for ODEModel
@test_approx_eq_eps evaluate!(m1,p1) datavalues(m1.measurements) 1.0
@test_approx_eq_eps evaluate!(m1,p1) springmassdata(t1,i1,[110.0,11.0]) 0.001
@test_approx_eq_eps evaluate!(m1,[110.0,11.0]) springmassdata(t1,i1,[110.0,11.0]) 0.001

###Test the loglikelihood function
@test_approx_eq loglikelihood(m1,r1) loglikelihood(m1.noisemodel,datavalues(m1.measurements),r1)

###Test the geometry function for base samples
s0 = samples(:base,2,1,Float64,Float64)
s0.values[1] = 100.0
s0.values[2] = 10.0
s1 = samples(:base,2,21,Float64,Float64)
s1.values[1,:] = 100.0
s1.values[2,:] = collect(9.5:0.05:10.5)
geometry!(m1,s0)
geometry!(m1,s1)

####################################################################################################

t2 = linspace(0.0,10.0,200)
d2 = [1.0,1.0]
p2 = [0.5,2.0/3.0]
v2 = [1e-2,1e-2]

m2 = sincosmodel(t2,p2,v2,d2)
r2 = evaluate!(m2,p2)

@test numparas(m2) == 2

###Test the evaluate function for TargetModel
@test_approx_eq_eps evaluate!(m2,p2) datavalues(m2.measurements) 2.0
@test_approx_eq_eps evaluate!(m2,p2) sincos(t2,p2) 0.001

@test_approx_eq loglikelihood(m2,r2) loglikelihood(m2.noisemodel,datavalues(m2.measurements),r2)

###Test the geometry function for base samples
s2 = samples(:base,2,21,Float64,Float64)
s2.values[1,:] = collect(0.49:0.001:0.51)
s2.values[2,:] = collect(2/3-0.01:0.001:2/3+0.0101)
geometry!(m2,s2)

####################################################################################################

t3 = linspace(0.0,10.0,200)
d3 = [1.0,1.0]
p3 = [0.5,2.0/3.0]
v3 = [1e-2,1e-2]

m3 = sincosmodel(t3,p3,v3,d3)
r3 = evaluate!(m3,p3)

@test numparas(m3) == 2

###Test the evaluate function for TargetModel
@test_approx_eq_eps evaluate!(m3,p3) datavalues(m3.measurements) 2.0
@test_approx_eq_eps evaluate!(m3,p3) sincos(t3,p3) 0.001

@test_approx_eq loglikelihood(m3,r3) loglikelihood(m3.noisemodel,datavalues(m3.measurements),r3)

###Test the geometry function for base samples
s3 = samples(:base,2,21,Float64,Float64)
s3.values[1,:] = collect(0.49:0.001:0.51)
s3.values[2,:] = collect(2/3-0.01:0.001:2/3+0.0101)
geometry!(m3,s3)

######################################################################################################

type UnitTestModelType <: AbstractModel end

@test_throws MethodError evaluate!(UnitTestModelType,zeros(3))

###test the show functions
println()
println("====================")
println("Test show() function")
println("====================")
show(m1)
show(m2)
println("====================")
println("End  show() function")
println("====================")
println()

















