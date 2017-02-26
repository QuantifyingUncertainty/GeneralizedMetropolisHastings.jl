# Set up a folder structure

## For a simple project

Prepare a separate `julia` package with the recommended folder structure:

```
MyGMHPackage.jl/
    MyGMHPackage.jl
    data/
    models/
    notebooks/
    scripts/
    test/
```

Measurement data sets are kept in `data`, model functions in `model`, and runnable MCMC
experiments in `scripts` and `notebooks`. Keep unit and other
tests for your package in `test`. Define other folders (e.g., for specific analysis code
or to store the results of the MCMC runs) as required.

`MyGMHPackage.jl` is the top level `module` file. It specifies imported packages,
included module files, and package exports. An example is available [here](https://github.com/QuantifyingUncertainty/GMHExamples.jl/blob/master/GMHExamples.jl).

## For a larger problem domain

When the problem domain of interest is likely to contain more than one type
of estimation experiment, repeat the above folder structure for different
experiment types, and group them hierarchically into problems and domains.

```
MyGMHPackage.jl/
    MyGMHPackage.jl
    domain1/
        problem1/
            data/
            models/
            notebooks/
            scripts/
        problem2/
            ...
        problem3/
            ...
    domain2/
        problem1/
            ...
        problem2/
            ...
    ...
    test/
```

Keep the `MyGMHPackage.jl` module file and `test` folder at the top level of the package.

An example of this folder structure can be seen
in [GMHExamples.jl](https://github.com/QuantifyingUncertainty/GMHExamples.jl.git).

## Add to Julia search path

Add the package to the `Julia` search path by adding the following line
to your `.juliarc.jl` file.

```julia
push!(LOAD_PATH,"/local/path/to/MyGMHPackage.jl")
```
