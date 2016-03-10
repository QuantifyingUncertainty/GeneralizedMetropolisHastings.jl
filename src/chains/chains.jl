abstract AbstractChain

chain{N<:Number,T<:AbstractFloat}(s::Symbol,nparas::Int,nsamples::Int,::Type{N},::Type{T},args...) = _chain(Val{s},nparas,nsamples,N,T,args...)

numparas(c::AbstractChain) = size(c.values,1)
numsamples(c::AbstractChain) = size(c.values,2)

accepted!(c::AbstractChain,indicator::AbstractIndicatorMatrix) = c.accepted += accepted(indicator)

samples(c::AbstractChain) = c.values
loglikelihood(c::AbstractChain) = c.loglikelihood
logprior(c::AbstractChain,m::AbstractModel) = logprior(parameters(m),c.values,eltype(c.loglikelihood))
logposterior(c::AbstractChain,m::AbstractModel) = c.loglikelihood + logprior(c,m)

