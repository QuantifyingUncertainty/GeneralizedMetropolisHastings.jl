i1 = indicator(:stationary,1,1,Float64)
i2 = indicator(:stationary,4,2,Float32)
i3 = indicator(:stationary,10,10,Float64)

@test numproposals(i1) == 1 && numsamples(i1) == 1 && length(i1.samples) == 2 && length(i1.stationary) == 2
@test numproposals(i2) == 4 && numsamples(i2) == 2 && length(i2.samples) == 3 && length(i2.stationary) == 5
@test numproposals(i3) == 10 && numsamples(i3) == 10 && length(i3.samples) == 11 && length(i3.stationary) == 11

u2 = log([0.5f0])
u3 = log([0.8])

v1 = log([0.6])
v1c = vcat(v1,[0.0])
v2 = Vector{Vector{Float32}}(1)
v2[1] = log([0.1f0,0.2f0,0.5f0,1.5f0])
v2c = vcat(v2[1]+u2[1],[0.0f0])
v3 = Vector{Vector{Float64}}(2)
v3[1] = log([1.0,1.5,2.0,0.5,0.3])
v3[2] = log([0.1,2.0,0.3,1.5,0.5])
v3c = vcat(v3[1]+u3[1],v3[2]+u3[1],[0.0])

a1 = transitionprobability!(i1,v1)
a2 = transitionprobability!(i2,u2,v2)
a3 = transitionprobability!(i3,u3,v3)

@test_approx_eq sum(a1) 1.0
@test_approx_eq a1 exp(v1c-maximum(v1c))./sum(exp(v1c-maximum(v1c)))
@test_approx_eq sum(a2) 1.0f0
@test_approx_eq a2 exp(v2c-maximum(v2c))./sum(exp(v2c-maximum(v2c)))
@test_approx_eq sum(a3) 1.0
@test_approx_eq a3 exp(v3c-maximum(v3c))./sum(exp(v3c-maximum(v3c)))

@test (srand(567) ; sampleindicator!(i1)) == (srand(567) ; vcat([2],[rand(Distributions.Categorical(a1)) for i=1:numsamples(i1)]))
@test (srand(567) ; sampleindicator!(i2)) == (srand(567) ; vcat([5],[rand(Distributions.Categorical(a2)) for i=1:numsamples(i2)]))
@test (srand(567) ; sampleindicator!(i3)) == (srand(567) ; vcat([11],[rand(Distributions.Categorical(a3)) for i=1:numsamples(i3)]))

@test accepted(i1) == sum(i1.samples[2:end].!=i1.samples[1:end-1])
@test accepted(i2) == sum(i2.samples[2:end].!=i2.samples[1:end-1])
@test accepted(i3) == sum(i3.samples[2:end].!=i3.samples[1:end-1])

###test the show functions
println()
println("====================")
println("Test show() function")
println("====================")
show(i1)
show(i2)
show(i3)
println("====================")
println("End  show() function")
println("====================")
println()
