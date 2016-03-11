## GeneralizedMetropolisHastings
Code base for the Generalized Metropolis-Hastings (GMH) algorithm for MCMC [(Calderhead, 2014)](#refs). 

Please check [GMH-Examples.jl](https://github.com/QuantifyingUncertainty/GMH-Examples.jl) for example models, scripts and notebooks using the GMH code base.

Instructions to install the package and run these experiments in Amazon Web Services, JuliaBox or on your local machine: http://quantifyinguncertainty.github.io

The **src** directory contains all the source code for the package.
	
The **test** directory contains tests of the package. You can run the tests from the Julia command line by executing:

```
julia> Pkg.test("GeneralizedMetropolisHastings")
```
###<a name="refs"/>References
Calderhead B. (2014), A general construction for parallelizing Metropolis-Hastings algorithms, PNAS, Vol: 111, Pages: 17408-17413
