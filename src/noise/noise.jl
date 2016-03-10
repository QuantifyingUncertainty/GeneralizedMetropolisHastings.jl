abstract AbstractNoiseModel

noise(s::Symbol,args...) = _noise(Val{s},args...)
noisemodelname(::AbstractNoiseModel) = "AbstractNoiseModel"

function show(io::IO,n::AbstractNoiseModel)
  println(io,noisemodelname(n))
  for f in fieldnames(n)
    print(io," ",f,": ",typeof(getfield(n,f)))
    println(io)
  end
  nothing
end

