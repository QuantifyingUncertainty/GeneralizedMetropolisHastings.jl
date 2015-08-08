import Distributions.Normal
import Distributions.Uniform
import GeneralizedMetropolisHastings.Prior
import GeneralizedMetropolisHastings.expand_names

#tests of the inner constructor
pa = ModelParameters{Uniform,ASCIIString}(["a","b"],fill(Uniform(0.0,10.0),2),[1.0,1.0])
pu = ModelParameters{Normal,ASCIIString}(fill("",2),[Normal(0.0,2.0),Normal(-2.0,5.0)],[1.0,1.0])
pn = ModelParameters{Uniform,UTF8String}(["μ1","μ2"],fill(Uniform(0.0,10.0),2),[1.0,1.0])
pp = ModelParameters{Prior,UTF8String}(["μ",""],Prior[Normal(0.0,1.0),Uniform(0.0,10.0)],[1.0,2.0])

@test pp.defaults[1] == pp.defaults[2]-1.0
@test pp.priors[1] != pp.priors[2]
@test named(pp) == [1] && anonymous(pp) == [2]

#tests of the function that inserts empty strings for unnamed parameters
@test expand_names(3,UTF8String[],Int[]) == fill(utf8(""),3)
@test expand_names(3,ASCIIString[],Int[]) == fill("",3)
@test expand_names(3,["a"],[2]) == ["","a",""]
@test expand_names(3,["μ","σ"],[1,3]) == ["μ","","σ"]
@test expand_names(3,["a","b","c"]) == ["a","b","c"]

#tests of the outer constructors
p1 = ModelParameters([1.0,2.0,3.0],[Uniform(0.0,1.0),Uniform(0.0,1.0),Uniform(0.0,1.0)];names=["a","μ"],index=[1,2])
p1a = ModelParameters([1.0,2.0,3.0],[Uniform(0.0,1.0),Uniform(0.0,1.0),Uniform(0.0,1.0)];names=["a","μ",""])
p1b = ModelParameters([1.0,2.0,3.0],[Uniform(0.0,1.0),Normal(0.0,1.0),Uniform(0.0,1.0)])

@test numel(p1) == 3 && p1.priors[1] == p1.priors[2] == p1.priors[3]
@test named(p1) == [1,2] && anonymous(p1) == [3]
@test named(p1b) == [] && anonymous(p1b) == [1,2,3]
@test typeof(p1b.priors[1]) == Uniform <: Prior
@test typeof(p1b.priors[2]) == Normal <: Prior

#tests of the equality operator
@test p1 == p1a && p1 != p1b && p1 != ModelParameters()

p2 = ModelParameters([1.0,2.0,3.0],Uniform(0.0,10.0);names=["a","b"],index=[1,2])
p2a = ModelParameters([1.0,2.0,3.0],Uniform(0.0,10.0);names=["a","b",""])
p2b = ModelParameters([0.0,0.0,0.0],Uniform(0.0,10.0))
p2c = ModelParameters([0.0,0.0,0.0])

@test p2 == p2a
@test p2c.priors[1] == Uniform(-1e43,1e43)

p3 = ModelParameters([Uniform(0.0,10.0),Uniform(0.0,10.0),Uniform(0.0,10.0)])
p4 = ModelParameters(3,Uniform(0.0,10.0))
p5 = ModelParameters(3)
p6 = ModelParameters(0)
p7 = ModelParameters()

@test p3 == p4
@test p2b == p4
@test p2c == p5
@test p2 != p4
@test p6 == p7
@test values(ValuesFromDefault(),p2) == [1.0,2.0,3.0]
@test values(ValuesFromDefault(),p3) == [0.0,0.0,0.0]
@test values(ValuesFromPrior(),p3) != values(ValuesFromPrior(),p4) #random values should not be the same
@test (srand(0);values(ValuesFromPrior(),p3)) == (srand(0);values(ValuesFromPrior(),p4)) #but they are the same if we reset the random generator

println("Test of the show() function")
println("====================")
println("Test show() function")
show(p1)
show(p1b)
show(p6)
println("End  show() function")
println("====================")
println()

nothing
