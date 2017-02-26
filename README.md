# GeneralizedMetropolisHastings

*A parallel Monte-Carlo Markov Chain package*

| **Documentation**   | **PackageEvaluator**                          |      **Build Status**        |
|:-------------------:|:---------------------------------------------:|:----------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url] | [![][pkg-0.4-img]][pkg-0.4-url] [![][pkg-0.5-img]][pkg-0.5-url] | [![][travis-img]][travis-url] [![][codecov-img]][codecov-url] |


Code base for the Generalized Metropolis-Hastings (GMH) algorithm [(Calderhead, 2014)](#refs).

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **most recently tagged version of the documentation.**
- [**LATEST**][docs-latest-url] &mdash; **in-development version of the documentation.**

For detailed instructions on how to set up and run GMH in Amazon Web Services or in JuliaBox see [Quantifying Uncertainty](http://quantifyinguncertainty.github.io)

Some elementary examples of how to use the code are provided in repository [GMHExamples.jl](https://github.com/QuantifyingUncertainty/GMHExamples.jl).

A biological model for use with the code is provided in repository [GMHPhotoReceptor.jl](https://github.com/QuantifyingUncertainty/GMHPhotoReceptor.jl)

###<a name="refs"/>References
Calderhead B. (2014), *A general construction for parallelizing Metropolis-Hastings algorithms*, PNAS, Vol: 111, Pages: 17408-17413 [10.1073/pnas.1408184111](http://www.pnas.org/content/111/49/17408.abstract)

[contrib-url]: https://quantifyinguncertainty.github.io/GeneralizedMetropolisHastings.jl/latest/man/contributing.html

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://quantifyinguncertainty.github.io/GeneralizedMetropolisHastings.jl/latest

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://quantifyinguncertainty.github.io/GeneralizedMetropolisHastings.jl/stable

[travis-img]: https://travis-ci.org/quantifyinguncertainty/GeneralizedMetropolisHastings.jl.svg?branch=master
[travis-url]: https://travis-ci.org/QuantifyingUncertainty/GeneralizedMetropolisHastings.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/egdu3hrptf3mnfc6/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/MichaelHatherly/GeneralizedMetropolisHastings-jl-bqgcw/branch/master

[codecov-img]: https://codecov.io/gh/quantifyinguncertainty/GeneralizedMetropolisHastings.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/quantifyinguncertainty/GeneralizedMetropolisHastings.jl

[issues-url]: https://github.com/quantifyinguncertainty/GeneralizedMetropolisHastings.jl/issues

[pkg-0.4-img]: http://pkg.julialang.org/badges/GeneralizedMetropolisHastings_0.4.svg
[pkg-0.4-url]: http://pkg.julialang.org/?pkg=GeneralizedMetropolisHastings
[pkg-0.5-img]: http://pkg.julialang.org/badges/GeneralizedMetropolisHastings_0.5.svg
[pkg-0.5-url]: http://pkg.julialang.org/?pkg=GeneralizedMetropolisHastings
[pkg-0.6-img]: http://pkg.julialang.org/badges/GeneralizedMetropolisHastings_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=GeneralizedMetropolisHastings
