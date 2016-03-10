abstract AbstractRunner

### Factory function
runner(s::Symbol,args...;keyargs...) = _runner(Val{s},args...;keyargs...)
