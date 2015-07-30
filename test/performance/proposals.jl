function testnormaldensity(m::Vector{Float64},c::Matrix{Float64},nprop::Int,imax::Int)

  println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  println("NormalDensity test with mean: ",m,", nprop=",nprop,", imax=",imax)

  #helper variables to test performance
  npars::Int = length(m)
  v1 = zeros(npars)
  v2 = zeros(npars,nprop)
  v3 = zeros(npars,nprop,imax)
  s = GeneralizedMetropolisHastings.BaseSample(npars)
  sp = [GeneralizedMetropolisHastings.BaseSample(npars) for i=1:nprop]

  #pre-allocate variables and arrays for testing
  n = Distributions.MvNormal(m,c)
  ni = Array(Distributions.MvNormal,imax)
  p = GeneralizedMetropolisHastings.NormalDensity(m,c)
  pp = [GeneralizedMetropolisHastings.NormalDensity(m,c) for i=1:nprop]
  pi= Array(GeneralizedMetropolisHastings.NormalDensity,imax)

  #performance test the constructors
  gc()
  println("===============================================")
  println("Construction of MvNormal from mean and covariance arrays")
  @time for i=1:imax
    ni[i] = Distributions.MvNormal(m,c)
  end
  gc()
  println("Construction of NormalDensity from mean and covariance arrays")
  @time for i=1:imax
    pi[i] = GeneralizedMetropolisHastings.NormalDensity(m,c)
  end

  #performance test of the update_density functions
  gc()
  println("================================================")
  println("Updating the mean of the density")
  @time for i=1:imax
    update_density!(pi[i],m)
  end
  println("Updating the covariance of the density")
  @time for i=1:imax
    update_density!(pi[i],c)
  end
  println("Updating the mean and the covariance of the density")
  @time for i=1:imax
    update_density!(pi[i],m,c)
  end

  #performance test the random sampling
  println("===============================================")
  println("Generating nprop*imax samples: (",nprop*imax,")")
  gc()
  println("Via ",imax," calls to rand! (vectorized version)")
  @time for i=1:imax
    rand!(n,v2)
  end
  gc()
  println("Via ",imax*nprop," calls to rand! (double for loop)")
  @time for i=1:imax
    for j=1:nprop
      rand!(n,v1)
    end
  end
  gc()
  println("===============================================")
  println("Via ",imax," calls to vectorized NormalDensity.propose!")
  @time for i=1:imax
    GeneralizedMetropolisHastings.propose!(p,sp)
  end
  gc()
  println("Via ",imax," calls to double vectorized NormalDensity.propose!")
  @time for i=1:imax
    GeneralizedMetropolisHastings.propose!(pp,sp)
  end
  println("===============================================")
  gc()

  nothing
end

#run the function to make sure it is compiled
testnormaldensity([4.0,3.0,2.0,1.0,0.0],eye(5),1,1)

nsamples = 100000
for nprop in [10,20,40]
  testnormaldensity([4.0,3.0,2.0,1.0,0.0],eye(5),nprop,div(nsamples,nprop))
end

nothing

