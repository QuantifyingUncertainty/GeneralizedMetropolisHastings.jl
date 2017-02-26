# GeneralizedMetropolisHastings

*A parallel Monte-Carlo Markov Chain package*

Code base for the Generalized Metropolis-Hastings (GMH) algorithm [(Calderhead, 2014)](#refs). 

The main documentation for the GMH package can be found here:
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://quantifyinguncertainty.github.io/GeneralizedMetropolisHastings.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://quantifyinguncertainty.github.io/GeneralizedMetropolisHastings.jl/latest)

For detailed instructions on how to set up and run MCMC experiments using GMH in Amazon Web Services, in JuliaBox or on your local machine, see: http://quantifyinguncertainty.github.io

Some elementary examples of how to use the code are provided in [GMHExamples.jl](https://github.com/QuantifyingUncertainty/GMHExamples.jl).

A biological model for use with the code is provided in [GMHPhotoReceptor.jl](https://github.com/QuantifyingUncertainty/GMHPhotoReceptor.jl)
###<a name="refs"/>References
Calderhead B. (2014), A general construction for parallelizing Metropolis-Hastings algorithms, PNAS, Vol: 111, Pages: 17408-17413
