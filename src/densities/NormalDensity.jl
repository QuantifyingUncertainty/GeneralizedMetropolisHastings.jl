type NormalDensity <: ProposalDensity
  normal::MvNormal
end

###Constructor with mean and covariance matrix
NormalDensity(μ::Vector{Float64},Σ::Matrix{Float64}) = NormalDensity(MvNormal(μ,Σ))
###Constructor with zero mean and covariance matrix
NormalDensity(Σ::Matrix{Float64}) = NormalDensity(MvNormal(zeros(size(Σ,1)),Σ))
###Constructor with zero mean and eye covariance matrix
NormalDensity(n::Int) = NormalDensity(zeros(n),eye(n))

###update the mean of the density
update_density!(d::NormalDensity,μ::Vector{Float64}) = @simd for i=1:length(d.normal.μ) @inbounds d.normal.μ[i] = μ[i] end

###update the covariance of the density
update_density!(d::NormalDensity,Σ::Matrix{Float64}) = (d.normal = MvNormal(mean(d.normal),Σ))

###update the mean and the covariance of the density
update_density!(d::NormalDensity,μ::Vector{Float64},Σ::Matrix{Float64}) = (d.normal = MvNormal(μ,Σ))

###draw a proposal and store it inside the values field of the sample
propose!(d::NormalDensity,s::MCSample) = rand!(d.normal,s.values)

###draw n proposals from the density
propose(d::NormalDensity) = rand(d.normal)
propose(d::NormalDensity,n::Int) = rand(d.normal,n)

###calculate the probability of a point, given a density
logprobability(d::NormalDensity,s::MCSample) = logpdf(d.normal,s.values)

###base functionality
==(d1::NormalDensity,d2::NormalDensity) = (mean(d1.normal) == mean(d2.normal) && cov(d1.normal) == cov(d2.normal))

function Base.show(io::IO,d::NormalDensity)
  println(io,"NormalDensity of length ",length(d.normal), " with")
  println(io,"mean: ",mean(d.normal))
  println(io,"covariance: ")
  println(io,cov(d.normal))
  println(io)
  nothing
end
