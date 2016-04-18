abstract AbstractRunner

### Factory function
runner(p::AbstractPolicy,args...;keyargs...) = _runner(p,args...;keyargs...)
