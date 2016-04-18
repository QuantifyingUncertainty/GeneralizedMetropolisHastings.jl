#Abstract runtime policy
abstract AbstractPolicy

#Factory function, currently supported type is :mh
policy(s::Symbol,args...;keyargs...) = _policy(Val{s},args...;keyargs...)
