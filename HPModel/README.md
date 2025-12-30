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
Here we leave some technical details of the functions and structures used in our code, it is advisable, to better understand the content of this document you are advised to read the pdf provided in the `docs/` directory. 
DISCLAIMER: The following functions may not be completely up to date, however the descriptions should be equivalent for most instances.

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
A `State` is nothing but the combination of a `Path` and an occupation matrix (`occupation`) this data structure was introduced to implement self avoidance in an efficient manner, however, memory allocation overhead, at least for short sequences turns out to be much more relevant than the cost of the slightly longer naive version. To avoid overloading the memory with potentially useless information we have decided to keep `State` and `Path` as separate structures that are used according to necessities.

### Core functions

#### Extention of a Self Avoiding Path
```julia
function extend_SAP(path::Path, dim::Int, base_vectors::Vector{Vector{Int}})::Vector{Path} 
    buffer = Vector{Path}(undef, 2*dim)
    last_position = path.vertices[path.length].position
    extensions_count = 0
    for axis_idx in 1:dim
        for sign in [+1,-1]
            step = sign * base_vectors[axis_idx]
            new_position = last_position + step
            if !is_in(path, new_position)
                extensions_count += 1
                buffer[extensions_count] = path_extend(path, new_position, dim)
            end
        end 
    end
    return buffer[1:extensions_count]
end
```
This function extends a path avoiding self intersection, in the source file you can find a function with the same name that extends `States` using the informtion contained in the occupation matrix

#### Equivalence of Paths
```julia
function are_equivalent(path1::Path, path2::Path, matrices)::Bool
    if path1.length == path2.length
        for matrix in matrices
            if are_equal(path1, linear_transform(path2, matrix))
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

    #1 
    if dim > 2
        max_size = Base.length(unique_SAPs(length,dim-1))*2
    else
        max_size = 1*2
    end
    #
    
    matrices = O(dim)
    seeds_buffer = Vector{State}(undef, max_size)
    states = Vector{State}(undef, max_size)
 
    next_states = Vector{State}(undef, max_size)
    next_states_candidates = Vector{State}(undef, max_size)

    states[1] = ZeroState(dim)
    states_count = 1

    seeds_count = 0
    next_states_count = 0

    #2
    for i in 1:length-1
        seeds[i] = seeds_buffer[1:seeds_count]
        next_states_count = 0
        next_states_candidates_count = 0
        seeds_count = 0

        #2.1
        for j in 1:states_count
            parent = states[j]
            children = extend_SAP(parent, dim)
            is_first_completely_asymmetric = true
            is_first_symmetry_break = true

            #2.1.1
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
        #    

        for j in 1:next_states_count
            states[j] = next_states[j]
        end
        states_count = next_states_count
    end
    #

    seeds_buffer[seeds_count+1:seeds_count+next_states_count] = next_states[1:next_states_count]
    seeds_count += next_states_count
    seeds[length] = seeds_buffer[1:seeds_count]
    return seeds
end
```
This function is the main implementation of the Seed-Sprout approach, it searches for seeds so that in the later stages they can be attached to sprouts. The details of the implementation are...

#### Search of Topological Paths

```julia
function topological_SAPs(depth::Integer, dim::Integer)::Vector{Vector{Path}}
    topological_SAPs = Vector{Vector{Path}}(undef, depth)
    topological_SAPs[1] = [ZeroPath(dim)]
    matrices = O(dim)
    base_vectors = compute_base_vectors(dim)
    for i in 1:depth-1
        buffer = [extend_SAP(path, dim, base_vectors) for path in topological_SAPs[i]]
        topological_SAPs[i+1] = fast_flatten([get_duplicate_free(vector, dim, matrices) for vector in buffer])
    end
    return topological_SAPs
end
```
#### Quotientation by any Linear Group

```julia
function get_duplicate_free(paths::Vector{Path}, dim::Integer, matrices)::Vector{Path}
    len = length(paths)
    symm_buffer = Vector{Path}(undef, len)
    asym_buffer = Vector{Path}(undef, len)
    symm_count = 0
    asym_count = 0
    for path in paths
        if path.asymmetry_flag & 1 == 1
            asym_count += 1
            asym_buffer[asym_count] = path
        else           
            is_unique = true
            for i in 1:symm_count
                is_unique &= !are_equivalent(symm_buffer[i], path, matrices)
            end
            if is_unique
                symm_count += 1
                symm_buffer[symm_count] = path
            end
        end
    end
    duplicate_free = Vector{Path}(undef, symm_count + asym_count)
    duplicate_free[1:symm_count] = symm_buffer[1:symm_count]
    duplicate_free[(symm_count+1):(symm_count+asym_count)] = asym_buffer[1:asym_count]
    return duplicate_free
end
```

#### Self Adjacency Triangle
```julia
function path_self_adjacency_triangle(path::Path)::Vector{Vector{Bool}}
    triangle = Vector([fill(false, i) for i in 1:path.length])
    for i in 3:path.length
        for j in 1:i-2
            if sum(abs.(path.vertices[i].position - path.vertices[j].position)) == 1
                triangle[i][j] = true
            end
        end
    end
    return triangle
end
```
#### Calculation of the self interaction (Any Model)
```julia
function calculate_self_interaction(sequence::Vector{Int}, triangle::Vector{Vector{Bool}}, interaction_model)
    interaction = 0
    for i in 1:Base.length(triangle)
        for j in 1:i
            interaction += triangle[i][j] * interaction_model(sequence[i], sequence[j])
        end
    end
    return interaction
end
```
#### Relevant Self Interations
```julia
function multiplicative_binary_interaction_model(a,b)
    return a*b
end

function lee_binary_interaction(a,b)
    if a == b == 1
        return 2.3
    elseif a == b == 0
        return 0
    else
        return 1
    end
end
```
#### Construction of all possible Binary sequence of a given length
```julia
function all_binary_sequences(len::Int)::Vector{Vector{Int}}
    sequences = Vector{Vector{Int}}([zeros(len) for _ in 1:2^len])
    for i in 1:2^len
        for j in 1:len
            sequences[i][j] = Int(UInt(i) >> (j-1) & UInt(1))
        end
    end
    return sequences
end
```
#### Construction of a random set of binary sequences
```julia
function random_binary_sequences(len::Int, number::Int)::Vector{Vector{Int}}
    sequences = Vector{Vector{Int}}([Vector{Int}(undef, len) for i in 1:number])
    for i in 1:number
        key = rand(1:2^len)
        for j in 1:len
            sequences[i][j] = Int(UInt(key) >> (j-1) & UInt(1))
        end
    end
    return sequences
end
```
#### Categorization by compactness
```julia
function categorize_by_compactness(paths::Vector{Path})::Vector{Vector{Path}}
    self_adjacency_triangles = path_self_adjacency_triangle.(paths)
    adjacency_count = count_adjacent.(self_adjacency_triangles)
    categories = Vector{Vector{Path}}([[] for _ in 1:maximum(adjacency_count)+1])
    for i in 1:Base.length(paths)
        push!(categories[adjacency_count[i]+1], path_copy(paths[i]))
    end
    return categories
end 
```
#### Calculation of g(t)
```julia
function calculate_g(sequences::Vector{Vector{Int}}, self_adjacency_triangles::Vector{Vector{Vector{Bool}}}, interaction_model)::Vector{Int}
    len = Base.length(self_adjacency_triangles)
    g = zeros(Int, len)
    for sequence in sequences
        Emax = 0
        multiplicity = 0
        for i in 1:len
            E = calculate_self_interaction(sequence, self_adjacency_triangles[i], interaction_model)
            if E == Emax
                multiplicity += 1
            elseif E > Emax
                multiplicity = 1
                Emax = E
            end
        end
    
        g[multiplicity] += 1
    end
    return g
end
```
### Additional Notes
Considering the scope, efficiency was prioretized over safety, many functions may throw unhandled errors.
