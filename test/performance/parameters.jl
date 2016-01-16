function create_parameters{T<:AbstractFloat}(s::Integer,::Type{T})
  randdefault = rand(T,s)
  randlow = rand(s) - 2
  randhigh = rand(s) + 2
  p1 = Vector{MCParameter}(s)
  p2 = Vector{MCParameterUnivariatePrior}(s)
  p3 = Vector{MCParameterDefaultValue}(s)
  gc() ; @time for i=1:2:s p1[i] = parameter(randdefault[s]) ; p1[i+1] = parameter(randlow[s],randhigh[s],randdefault[s]) end
  gc() ; @time for i=1:s p2[i] = parameter(randlow[s],randhigh[s],randdefault[s]) end
  gc() ; @time for i=1:s p3[i] = parameter(randdefault[s]) end
  println("Memory taken by ",typeof(p1[1]),": ",objectsizetostr(p1))
  println("Memory taken by ",typeof(p2[1]),": ",objectsizetostr(p2))
  println("Memory taken by ",typeof(p3[1]),": ",objectsizetostr(p3))
end

create_parameters(100000,Float64)
create_parameters(100000,Float32)
