abstract AbstractChain

chain{N<:Number,T<:AbstractFloat}(s::Symbol,nparas::Int,nsamples::Int,::Type{N},::Type{T},args...;keyargs...) = _chain(Val{s},nparas,nsamples,N,T,args...;keyargs...)

numparas(c::AbstractChain) = size(c.values,1)
numsamples(c::AbstractChain;store=:all) = _numsamples(Val{store},c)
numaccepted(c::AbstractChain) = c.numaccepted
numproposed(c::AbstractChain) = c.numproposed

accepted!(c::AbstractChain,indicator::AbstractIndicatorMatrix) = c.numaccepted += accepted(indicator)

samples(c::AbstractChain;paras=1:numparas(c),store=:all) = view(c.values,paras,_index(c,store))
loglikelihood(c::AbstractChain;store=:all) = view(c.loglikelihood,_index(c,store))
logprior(c::AbstractChain,m::AbstractModel;store=:all) = logprior(parameters(m),samples(c;store=store),eltype(c.loglikelihood))
logposterior(c::AbstractChain,m::AbstractModel;store=:all) = loglikelihood(c;store=store) + logprior(c,m;store=store)

mean(c::AbstractChain;store=:all) = mean(samples(c;store=store),2)
median(c::AbstractChain;store=:all) = numsamples(c;store=store)>0?median(samples(c;store=store),2):NaN
std(c::AbstractChain;store=:all) = std(samples(c;store=store),2)
var(c::AbstractChain;store=:all) = var(samples(c;store=store),2)

###Helper functions to index the results
@inline _index(c::AbstractChain,store::Symbol) = _startindex(Val{store},c):_endindex(Val{store},c)

@inline _startindex(::Type{Val{:burnin}},c::AbstractChain) = 1
@inline _startindex(::Type{Val{:main}},c::AbstractChain) = c.numburnin+1
@inline _startindex(::Type{Val{:proposed}},c::AbstractChain) = 1
@inline _startindex(::Type{Val{:all}},c::AbstractChain) = 1

@inline _endindex(::Type{Val{:burnin}},c::AbstractChain) = _numsamples(Val{:burnin},c)
@inline _endindex(::Type{Val{:main}},c::AbstractChain) = _numsamples(Val{:burnin},c) + _numsamples(Val{:main},c)
@inline _endindex(::Type{Val{:proposed}},c::AbstractChain) = _numsamples(Val{:proposed},c)
@inline _endindex(::Type{Val{:all}},c::AbstractChain) = _numsamples(Val{:all},c)

@inline _numsamples(::Type{Val{:burnin}},c::AbstractChain) = min(c.numproposed,c.numburnin)
@inline _numsamples(::Type{Val{:main}},c::AbstractChain) = max(0,c.numproposed-c.numburnin)
@inline _numsamples(::Type{Val{:proposed}},c::AbstractChain) = c.numproposed
@inline _numsamples(::Type{Val{:all}},c::AbstractChain) = length(c.loglikelihood)
