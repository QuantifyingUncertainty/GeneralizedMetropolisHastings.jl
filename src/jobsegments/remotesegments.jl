abstract AbstractRemoteSegments

remotesegments(p::AbstractPolicy,m::AbstractModel,s::AbstractSampler,nproposals::Int) = _remotesegments(p,m,s,nproposals)

numsegments(s::AbstractRemoteSegments) = s.numsegments
numproposalspersegment(s::AbstractRemoteSegments) = s.numproposalspersegment
numtotalproposals(s::AbstractRemoteSegments) = s.numsegments*s.numproposalspersegment

@inline _numjobsegments(::Type{Val{:procs}}) = nprocs()
@inline _numjobsegments(::Type{Val{:workers}}) = nworkers()
@inline _numjobsegments(::Type{Val{:none}}) = 0
@inline _numjobsegments(::Type{Val{:one}}) = 1
@inline _numjobsegments(::Type{Val{:two}}) = 2

@inline _processnumbers(::Type{Val{:procs}}) = procs()
@inline _processnumbers(::Type{Val{:workers}}) = workers()
@inline _processnumbers(::Type{Val{:none}}) = []
@inline _processnumbers(::Type{Val{:one}}) = workers()
@inline _processnumbers(::Type{Val{:two}}) = workers()

@inline _numproposalspersegment(nproposals::Int,njobsegments::Int) = ceil(Int,nproposals/njobsegments)
