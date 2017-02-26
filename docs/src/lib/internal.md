# Internal Documentation

```@meta
CurrentModule = GeneralizedMetropolisHastings
```

Documentation for `GeneralizedMetropolisHastings.jl`'s internal functions.

You may need to extend some of these types or functions when developing new
problem domains or extending the GMH package.

See [Public Documentation](@ref) for its public interface.

## Contents

```@contents
Pages = ["internal.md"]
Depth = 3
```

## Index

```@index
Pages = ["internal.md"]
```

## Internal functions and types

### Data - Internal

#### Data - Internal Functions
```@docs
datatypename(::AbstractData)
eltype(::AbstractData)
generate!(::AbstractData)
```

#### Data - Internal Types
```@docs
DataArray
DataFunction
```
