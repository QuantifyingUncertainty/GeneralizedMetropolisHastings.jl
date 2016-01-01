type NormalDensity end
#type SymmetricProposal <: ProposalDensity
#  density::Distributions
#end

###draw a proposal and store it inside the values field of the sample
#propose!(d::SymmetricProposal,s::MCSample) = (rand!(d.density,s.values) ; s)

#acceptance!(d::SymmetricProposal,s::MCSample) =


# function Base.show(io::IO,d::NormalDensity)
#   println(io,"NormalDensity of length ",length(d.normal), " with")
#   println(io,"mean: ",mean(d.normal))
#   println(io,"covariance: ")
#   println(io,cov(d.normal))
#   println(io)
#   nothing
# end
