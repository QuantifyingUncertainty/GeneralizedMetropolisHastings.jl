testpdf1!(r::AbstractVector,d::Distribution,m::Matrix{Float64}) = (logpdf!(r,d,m) ; r)
testpdf2!(r::AbstractVector,d::Distribution,m::Matrix{Float64}) = (for i=1:length(r) r[i] = logpdf(d,m[:,i]) end ; r)
testpdf3!(r::AbstractVector,d::Distribution,v::Vector{Vector{Float64}}) = (for i=1:length(r) r[i] = logpdf(d,v[i]) end ; r)
