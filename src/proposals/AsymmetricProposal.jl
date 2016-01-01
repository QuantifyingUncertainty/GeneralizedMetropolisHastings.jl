# type LogNormalDensity <: AsymmetricProposalDensity
#   μ::Vector{Float64} #location of the log-normal distribution
#   Σ::Matrix{Float64} #scale of the log-normal distribution
#   normal::Distributions.MvNormal #standard normal variate
#   LogNormalDensity(μ::Vector{Float64},Σ::Matrix{Float64}) = new(μ,Σ,Distributions.MvNormal(zeros(μ),1.0))
# end

# ###Constructor with zero mean and covariance matrix
# LogNormalDensity(Σ::Matrix{Float64}) = LogNormalDensity(zeros(size(Σ,1)),Σ)
# ###Constructor with zero mean and eye covariance matrix
# LogNormalDensity(n::Int) = LogNormalDensity(zeros(n),eye(n))

# ###update the mean of the density
# update_density!(d::LogNormalDensity,μ::Vector{Float64}) = (copy!(d.μ,μ) ; d)

# ###update the covariance of the density
# update_density!(d::LogNormalDensity,Σ::Matrix{Float64}) = (copy!(d.Σ,Σ) ; d)

# ###update the mean and the covariance of the density
# update_density!(d::LogNormalDensity,μ::Vector{Float64},Σ::Matrix{Float64}) = (copy!(d.μ,μ) ; copy!(d.Σ,Σ) ; d)

# ###draw a proposal and store it inside the values field of the sample
# propose!(d::LogNormalDensity,s::MCSample) = (rand!(d.normal,s.values) ; s.values = exp(d.μ + d.Σ*s.values) ; s)

# ###calculate the probability of a point, given a density
# logprobability(d::LogNormalDensity,s::MCSample) = logpdf(d.normal,\(d.Σ,log(s.values)-d.μ))

# ###base functionality
# import Base.==
# ==(d1::LogNormalDensity,d2::LogNormalDensity) = (d1.μ == d2.μ && d1.Σ == d2.Σ)

# function Base.show(io::IO,d::LogNormalDensity)
#   println(io,"LogNormalDensity of length ",length(d.μ), " with")
#   println(io,"location: ",d.μ)
#   println(io,"scale: ")
#   println(io,d.Σ)
#   println(io)
#   nothing
# end
