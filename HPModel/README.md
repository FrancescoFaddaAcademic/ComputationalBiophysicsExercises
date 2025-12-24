# HPModel.jl
HPModel is a Julia pakage that provides tools for the 
## Installation
To install the package in your own project run in the repl
```julia
  using Pkg
  Pkg.add(url=)
```
## Documentation
```julia

```
### Core structures
```julia
mutable struct Path{V<:AbstractVertex}<:AbstractPath
    vertices::Vector{V}
    length::Int
    asymmetry_flag::UInt
end
```
`Path` is synonimous to conformation in biological terms, a path contains an array `vertices` of vertices (that depending on their concrete type may contain different informations), their `length` (for efficient data transfer) and an `asymmetry_flag` that contains useful information on the degree of asymmetry of a certain path.  

```julia
mutable struct State
    path::Path
    occupation::Array{Bool}
end
```
A `State` is nothing but the combination of a `Path` and an occupation matrix (`occupation`) this data structure is useful to implement self avoidance in an efficient manner, to avoid overloading the memory with potentially useless information we have decided to keep `State` and `Path` as separate structures that are used according to necessities.

### Core functions

```julia

```
### Architecture
