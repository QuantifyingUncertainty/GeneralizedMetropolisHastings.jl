function testmetropolis(m::Vector{Float64},c::Matrix{Float64},nprop::Int,imax::Int)

  println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  println("MetropolisHastings test with normal density: ",m,", nprop=",nprop,", imax=",imax)

  #helper variables to test performance
  npars::Int = length(m)
  n = GeneralizedMetropolisHastings.NormalDensity(m,c)
  s = GeneralizedMetropolisHastings.MHNormal(npars,1.0)
  f = GeneralizedMetropolisHastings.BaseSample(m)

  #pre-allocate variables and arrays for testing
  b = [GeneralizedMetropolisHastings.BaseSample(npars) for i=1:nprop]
  h = GeneralizedMetropolisHastings.MHHeap(s,nprop)

  #performance test setting of the from field
  gc()
  println("===============================================")
  println("Via ",imax," calls to from(sample)!")
  @time for i=1:imax
    GeneralizedMetropolisHastings.set_from!(s,h,f)
  end
  gc()
  #performance test the random sampling, comparing it to sampling straight from NormalDensity
  gc()
  println("===============================================")
  println("Via ",imax," calls to vectorized NormalDensity.propose!")
  @time for i=1:imax
    GeneralizedMetropolisHastings.propose!(n,b)
  end
  gc()
  println("Via ",imax," calls to MMHeap.propose!")
  @time for i=1:imax
    GeneralizedMetropolisHastings.propose!(s,h,1,f)
  end
  println("===============================================")
  gc()

  nothing
end

#run the function to make sure it is compiled
testmetropolis([4.0,3.0,2.0,1.0,0.0],eye(5),2,1)

nsamples = 100000
for nprop in [10,20,40]
  testmetropolis([4.0,3.0,2.0,1.0,0.0],eye(5),nprop,div(nsamples,nprop))
end

nothing

