###test a Data object with one variable, and generating Float32 values
d0 = data(:array,0.0f0:0.1f0:1.0f0,zeros(Float32,11))
@test numvalues(d0) == 11
@test numvars(d0) == 1
@test generate!(d0).values == zeros(Float32,11)
@test dataindex(d0) == collect(0.0f0:0.1f0:1.0f0)
@test datavalues(d0) == zeros(Float32,11)
@test eltype(d0) == Float32

###measurment time
t1 = 0.0:0.1:1.0

###data generating functions
f1(t,a) = hcat(sin(a*t),cos(a*t))
f1!(r,t,a) = (for i=1:length(t) r[i,:] = [sin(a*t[i]),cos(a*t[i])] end ; r)
f2!(r,μ::Float64,Σ::Float64) = (d = Distributions.Normal(μ,Σ) ; rand!(d,r) ; r)

###data array
a1 = hcat(sin(2*pi*t1),cos(2*pi*t1))

###create the Data objects
d1 = data(:array,t1,a1)
d2 = data(:function,t1,f1,t1,2π)
d2! = data(:function!,t1,zeros(length(t1),2),f1!,t1,2π)
d3! = data(:function!,t1,zeros(length(t1),3),f2!,1.0,0.1)
d4! = data(:function!,t1,zeros(length(t1),3),rand!)

@test numvalues(d1) == numvalues(d2) == numvalues(d2!) == numvalues(d3!) == numvalues(d4!) == 11
@test numvars(d1) == numvars(d2) == numvars(d2!) == 2 && numvars(d3!) == numvars(d4!) == 3
@test datavalues(d1) == datavalues(d2) == datavalues(d2!) == a1
@test datavalues(generate!(d1)) == datavalues(generate!(d2)) == datavalues(generate!(d2!)) == a1
@test dataindex(d1) == dataindex(d2) == dataindex(d2!) == dataindex(d3!) == dataindex(d4!) == collect(t1)
@test eltype(d1) == eltype(d2) == eltype(d2!) == eltype(d3!) == eltype(d4!) == Float64
@test (srand(5643) ; datavalues(generate!(d3!))) == (srand(5643) ; rand(Distributions.Normal(1.0,0.1),11,3)) == datavalues(d3!)
@test (srand(5643) ; datavalues(generate!(d4!))) == (srand(5643) ; rand(11,3)) == datavalues(d4!)

###test obtaining a new copy of the datavector
e1 = datavalues(d1) ; e1[1] = 10.0
e2 = datavalues(d2) ; e2[1] = 10.0
@test e1 == e2 && e1 == datavalues(d1) && e2 == datavalues(d2)
@test e1 == datavalues(generate!(d1)) && e2 == datavalues(generate!(d2))

###noise model
s1 = [0.1,0.2]
n1 = noise(:gaussian,s1)

###apply a noise model to model data
c1 = datavalues(d1)
c2 = datavalues(d2)
srand(48032) ; m1 = copy(datavalues(d1)) ; for i=1:numvars(d1) m1[:,i] += rand(n1.distributions[i],numvalues(d1)) end ; m1
srand(48032) ; m2 = applynoise!(n1,copy(datavalues(d1)))
@test_approx_eq c1 datavalues(d1)
@test_approx_eq c2 datavalues(d2)
@test_approx_eq m1 m2

###test the loglikelihood of the data relative to the noise model
@test_approx_eq loglikelihood(n1,m1,datavalues(d1)) sum(broadcast((μ,σ,x)->Distributions.logpdf(Distributions.Normal(μ,σ),x),datavalues(d1),sqrt(s1)',m1))

println()
println("====================")
println("Test show() function")
show(d1)
show(d2)
show(d3)
show(d4)
show(n1)
println("End  show() function")
println("====================")
println()

nothing
