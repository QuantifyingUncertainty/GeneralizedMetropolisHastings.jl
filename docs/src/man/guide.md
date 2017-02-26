```@meta
CurrentModule = GeneralizedMetropolisHastings
```
# Package Guide

## Installation Guidelines

### Installing the GMH package to replicate published experiments

To run the GMH package with published experiments and data sets on
Amazon Web Services or JuliaBox, see the detailed [installation instructions](https://QuantifyingUncertainty.github.io)

### Installing the GMH package locally

The package is registered in `METADATA.jl` and can be installed using `Pkg.add`.

```julia
Pkg.add("GeneralizedMetropolisHastings")
```

You can test if the package is working correctly by running

```julia
Pkg.test("GeneralizedMetropolisHastings")
```

## Running existing experiments

A set of ready-to-run, elementary examples is available in [GMHExamples.jl](https://github.com/QuantifyingUncertainty/GMHExamples.jl.git).

Information on how to replicate published experiments in specific problem domains
is available from [Quantifying Uncertainty](https://QuantifyingUncertainty.github.io).

## Writing new experiments in existing problem domains

TBD

## Developing new problem domains

The following steps need to be performed to specify a new problem domain.

* [Set up a folder structure](@ref)
* [Add measurement data](@ref)
* [Specify measurement noise](@ref)
* [Specify model functions](@ref)

If the above steps have been completed, runnable MCMCM experiments
can be defined in `julia` scripts or `IJulia` notebooks as outlined in
[Writing new experiments in existing problem domains](@ref).
