# HPModel.jl
HPModel.jl is a Julia pakage that provides tools and examples relative to the Hydrophobic Polar Model. 
The directory `test/` contains two scripts that reproduce the required results of the two papers `Chan-Dill-1989` and `Lau-Dill-1989`  

### Implemented extensions
The pakage supports:
- Any dimension
- Any first neighbor interaction model that is compatible with the implemented data structures. In particular it supports the interaction model proposed by lee et al.
- Computation of random sequences for lengthy structures
(The compactness is calculated in the plot from `Chan-Dill-1989`)

## Documentation
Here we leave some technical details of the functions and structures used in our code, it is advisable, to better understand the content of this document, to first read the [].
### Core structures
#### Path
```julia
mutable struct Path{V<:AbstractVertex}<:AbstractPath
    vertices::Vector{V}
    length::Int
    asymmetry_flag::UInt
end
```
`Path` is synonimous to conformation in biological terms, a path contains an array `vertices` of vertices (that depending on their concrete type may contain different informations), their `length` (for efficient data transfer) and an `asymmetry_flag` that contains useful information on the degree of asymmetry of a certain path.  
#### State
```julia
mutable struct State
    path::Path
    occupation::Array{Bool}
end
```
A `State` is nothing but the combination of a `Path` and an occupation matrix (`occupation`) this data structure is useful to implement self avoidance in an efficient manner, to avoid overloading the memory with potentially useless information we have decided to keep `State` and `Path` as separate structures that are used according to necessities.

### Core functions

#### Extention of a Self Avoiding Path
```julia
function extend_SAP(state::State, dim::Int)::Vector{State} 
    path = state.path
    occupation = state.occupation
    new_state_buffer = Vector{State}(undef, 2*dim)
    last_position = path.vertices[path.length].position
    buffer_idx = 1
    extensions_count = 0
    base = compute_base_vectors(dim)
    for axis_idx in 1:dim
        for sign in [+1,-1]
            step = sign * base[axis_idx]
            new_position = last_position + step
            centered_idx = centered_index(new_position, path.length)
            if !occupation[centered_idx...]
                new_state_buffer[buffer_idx] = state_extend(state, new_position, dim)
                buffer_idx += 1
                extensions_count += 1
            end
        end 
    end
    return new_state_buffer[1:extensions_count]
end
```
This function extends a path avoiding self intersection, here the `State` structure comes useful to efficiently check the occupation of positions that are adjacent to that of the seed

#### Equivalence of Paths
```julia
function are_equivalent(path1::Path, path2::Path, matrices)::Bool
    if path1.length == path2.length
        for matrix in matrices
            if are_equal(long_path, linear_transform(path1, matrix))
                return true
            end
        end
        return false
    end
    return false
end
```
This function checks if two paths can be superimposed by virtue of one of transformations provided in the argument as matrices, in our case these are nothing but the orthogonal matrices on $\mathbb Z_N$.

#### Seed search
```julia
function search_seeds(length::Int, dim::Int)::Vector{Vector{State}}          
    seeds = Vector{Vector{State}}(undef, length) 
    if dim > 2
        max_size = Base.length(unique_SAPs(length,dim-1))*2
    else
        max_size = 1*2
    end
    matrices = O(dim)
    seeds_buffer = Vector{State}(undef, max_size)
    states = Vector{State}(undef, max_size)
 
    next_states = Vector{State}(undef, max_size)
    next_states_candidates = Vector{State}(undef, max_size)

    states[1] = ZeroState(dim)
    states_count = 1

    seeds_count = 0
    next_states_count = 0
    for i in 1:length-1
        seeds[i] = seeds_buffer[1:seeds_count]
        next_states_count = 0
        next_states_candidates_count = 0
        seeds_count = 0
        for j in 1:states_count
            parent = states[j]
            children = extend_SAP(parent, dim)
            is_first_completely_asymmetric = true
            is_first_symmetry_break = true
            for child in children
                ica = is_completely_asymmetric(child.path, dim)
                symmetry_break = child.path.asymmetry_flag != parent.path.asymmetry_flag
                if ica && is_first_completely_asymmetric
                    seeds_count += 1
                    seeds_buffer[seeds_count] = state_copy(child)
                    is_first_completely_asymmetric = false
                elseif symmetry_break && is_first_symmetry_break && !ica
                    next_states_count += 1
                    next_states[next_states_count] = state_copy(child)
                    is_first_symmetry_break = false
                elseif !((symmetry_break && !is_first_symmetry_break) || (ica && !is_first_completely_asymmetric))
                    next_states_candidates_count += 1
                    next_states_candidates[next_states_candidates_count] = state_copy(child)                    
                end
            end
        end

        if next_states_candidates_count != 0
            new_next_states = filter_equivalent(next_states_candidates[1:next_states_candidates_count], dim, matrices)
            new_next_states_count = Base.length(new_next_states)
            next_states[next_states_count+1:next_states_count + new_next_states_count] = new_next_states
            next_states_count += new_next_states_count
        end

        for j in 1:next_states_count
            states[j] = next_states[j]
        end
        states_count = next_states_count
    end
    seeds_buffer[seeds_count+1:seeds_count+next_states_count] = next_states[1:next_states_count]
    seeds_count += next_states_count
    seeds[length] = seeds_buffer[1:seeds_count]
    return seeds
end
```


### Additional Notes
Considering the scope, efficiency was prioretized over safety, many functions may throw unhandled errors.
