## GeneralizedMetropolisHastings
Code base for the Generalized Metropolis-Hastings (GMH) algorithm for MCMC [(Calderhead, 2014)](#refs). 

For instructions on how to set up and run MCMC experiments using GMH in Amazon Web Services, in JuliaBox or on your local machine, see: http://quantifyinguncertainty.github.io

Use the code with the models provided in [GMHModels.jl](https://github.com/QuantifyingUncertainty/GMHModels.jl) and the example scripts and notebooks provided in [GMHExamples.jl](https://github.com/QuantifyingUncertainty/GMH-Examples.jl)

From inside Julia, you can easily install the package using the following command:

```
julia> Pkg.clone("git://github.com/QuantifyingUncertainty/GeneralizedMetropolisHastings.jl")
```

This will install the package into your Julia package directory. 

The **src** directory contains all the source code for the package.
	
The **test** directory contains tests of the package. You can run the tests from the Julia command line by executing:

```
julia> Pkg.test("GeneralizedMetropolisHastings")
```
###<a name="refs"/>References
Calderhead B. (2014), A general construction for parallelizing Metropolis-Hastings algorithms, PNAS, Vol: 111, Pages: 17408-17413
