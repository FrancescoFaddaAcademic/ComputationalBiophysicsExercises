# HPModel

## Brief Theoretical Introduction and Architecture of the code
The Hydrophobic-Polar (HP) model is a toy model used to understand how the polar characteristics of a protein may affect its folding and conformation. In particular it is interesting in showing how the compactness of a protein is associated to a small number of ground state configurations 

### Implemented extensions
The pakage supports:
- Any dimension
- Any first neighbor interaction model that is compatible with the implemented data structures. In particular it supports the interaction model proposed by lee et al.
- Computation of random sequences for lengthy structures
(The compactness is calculated in the plot from `Chan-Dill-1989`)

### The Mathematics of Self Avoiding Walks
The study of topological configuration of self avoiding walks (or paths) on any n-dimensional grid is a well established area of mathematical research, for this reason, most of the nomenclature will be inspired to the corresponding mathematical concepts. In particular I will refer to Self Avoding Paths (SAP) instead of Secondary Structures, in the biological setting, these can be taken as synonims.
Let's start by considering the set of all possible (non self avoiding) paths pinned at the origin of an D-dimensional grid. It is easy to calculate the number of elements for each lentgh (as there are no collisions) as $(2D)^L$ (here the length here does not take into account the first point at the origin, however in the code it will be considered) This set has more structure to it, in particular there is a very natural graph structure that arises by thinking about each configuration as a node, and by letting two nodes be connected iff they differ by a last step. Moreover this graph structure is a tree (it has no loops). The reason behind this is that if we branch from a common configuration there is no possible way to reconnect later, as modifications can be done only at the end of the path. 

The tree of self avoiding pahts is nothing but the subtree where all self colliding branches have been pruned. We are not quite there yet, we would like to ignore orientation and (in some case mirror symmetry) the general idea behind this is the quotienting of a set under the action of a group, in our particular case the group will be the rigid linear symmetry group of the D-square-lattice and our set will be the tree of SAPs. Two SAPs are equivalent if they are connected by the action of some element of the group (if they are part of the same orbit). The set of orbits is what we are interested in for our understanding of topological SAPs. This quotient space is still a tree, however this is a little harder to see. In the following illustration I will give an intuitive idea of how this quotienting can be thougt of as a sort of folding of the original graph 

The natural approach for the computation of all the possible topological configurations is to computate of all possible SAPs and then check if they are in the same orbit of another element that we have already counted. However (while not in use) i've implemented another approach that looks at the roots of those subgraphs that contan SAPs that do not contain equivalent elements. I've named these roots Seeds, they are represented in light blue in the graph that I've drawn. The interesting property of these seeds is that, no matter what is attached to them, the composite SAP will always be asymmetric. In a sense they are minimal generators of complete asymmetry. What this means computationally is that it is sufficient to calculate them uniquely (removing all their copies) and then attach all possible SAPs that complete the length of the seed to that required for the calculation, we've called these completions sprouts. It is easy to show that the smallest seed for any dimension must have the dimension as its length. What this implies is that the biggest sprout one needs to calculate is the final length minus the dimensionality of the system. This may seem small but due too the exponential nature of the problem may have a big impact on the depth of the computation, mostly for the higher dimensional case. An objection that needs to be addressed is the following: The computation of these seeds may be computationally expensive, what tells us that this approach might be efficient? Here, while I do not have a definite answer, i would like to note that the seeds grow in number as the number of (D-1) SAPs of length lower than that of the problem in consideration. in particular for dimension 2 they grow linearly with length as for each length in 1 dimension there is a unique topological SAP. Again due to the exponential nature of the computation this contribution will always be vanishing with respect to the computation of the sprouts.
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
This function extends a path avoiding self intersection, in particular it tries to add to the input path all possible extensions by attaching a vertex, accepts or reject the configuration by checking a collision. In the source file you can find a function with the same name that extends `States` using the informtion contained in the occupation matrix.

#### Equivalence of Paths
```julia
function are_equivalent(path1::Path, path2::Path, matrices)::Bool
    if path1.length == path2.length
        for matrix in matrices
            linear_transform!(path2, matrix)
            if are_equal(path1, path2)
                linear_transform!(path2, inv(matrix))
                return true
            end
            linear_transform!(path2, inv(matrix))
        end
    end
    return false
end
```
This function checks if two paths can be superimposed by virtue of one of transformations provided in the argument as matrices, in our case these are nothing but the orthogonal matrices on $\mathbb Z_N$. In particular the method takes two paths, and after checking that their length is equal, transforms the second, compares it to the first one and, importantly, transforms it back. We must note that the `linear_transform!` function modifies the content of the path without creating a copy (this was done to avoid unnecessary slow memory alloactions), this implies that one needs to put everything back to where it was in order to avoid strange behaviour in mutable structs that may reference the same memory.

#### Quotientation by any Group (representation)

```julia
function get_duplicate_free(paths::Vector{Path}, dim::Integer, matrices)::Vector{Path}
    len = length(paths)
    buffer = Vector{Path}(undef, len)
    buffer_count = 0
    for path in paths      
        is_unique = true
        for i in 1:buffer_count
            is_unique &= !are_equivalent(buffer[i], path, matrices)
        end
        if is_unique
            buffer_count += 1
            buffer[buffer_count] = path
        end
    end
    duplicate_free = Vector{Path}(undef, buffer_count)
    duplicate_free[1:buffer_count] = buffer[1:buffer_count]
    return duplicate_free
end
```
This functions filters a set of paths in such a way that only one representative of any equivalence class is preserved, it uses the `are_equivalent` function to check if the current path is topologically equivalent to another, already saved instance.

#### Search of Topological Paths

```julia
function topological_SAPs(depth::Integer, dim::Integer)::Vector{Vector{Path}}
    topological_SAPs = Vector{Vector{Path}}(undef, depth)
    topological_SAPs[1] = [ZeroPath(dim)]
    matrices = O(dim)
    base_vectors = compute_base_vectors(dim)
    for i in 1:depth-1
        buffer = Vector{Path}(undef, (2dim)*length(topological_SAPs[i]))
        buffer_count = 0
        for path in topological_SAPs[i]
            extended_SAPs = extend_SAP(path, dim, base_vectors)
            if (path.asymmetry_flag ⊻ asymm_key(dim)) >> 1 == 0
                len = length(extended_SAPs)
                buffer[buffer_count+1:buffer_count + len] = extended_SAPs
            else 
                duplicate_free = get_duplicate_free(extended_SAPs, dim, matrices)
                len = length(duplicate_free)
                buffer[buffer_count+1:buffer_count + len] = duplicate_free
            end
            buffer_count += len
        end
        topological_SAPs[i+1] = buffer[1:buffer_count]
    end
    return topological_SAPs
end
```
This is the main outer implementation of the topological paths search. It constructs iteratively the set of self avoiding paths of length lesser or equal to `depth`. It starts from the previous order, produces all the extensions, checks the asymmetry of the parent path, and, in case this is not maximal, filters duplicates using the previous function, then, finally stores the elements.



#### Self Adjacency Triangle
```julia
function path_self_adjacency_triangle(path::Path)::Vector{Vector{Bool}}
    len = path.length
    triangle = Vector{Vector{Bool}}(undef, len)
    @inbounds for i in 1:len
        triangle[i] = Vector{Bool}(undef, i)
        @inbounds for j in 1:i
            triangle[i][j] = false
        end
    end
    dim = length(path.vertices[1].position)
    for i in 3:path.length
        for j in 1:i-2
            if one_distance(path.vertices[i], path.vertices[j], dim) == 1
                triangle[i][j] = true
            end
        end
    end
    return triangle
end
```
This funciton computes the self adjacency triangle of a given path, the structure has been chosen as a triangle instead of a matrix to reduce memory usage, however, for better performances (locality of data) it would be advisable to use a traditional matrix.

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
This function uses the previously calculated adjacency triangle, together with an arbitrary interaction model and calculates the cumulative self interaction of the entire path

#### Relevant Self Interactions
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
For example these are the two interactions that I've included in the code, they both are meant to work on binary sequences, however, any generalizaztion to a non-binary model should work correctly. 

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
````
This algorithm calculates all the possible binary sequences using some very efficient bitwise operations. In particular one can think of any number as a binary sequence, by iterating between 1 and 2^len one can easily produce all possible binary sequence, then it is just a matter of selecting the correct bit (by shifting it on the rightside position and masking it with 1) and inserting the value in the corresponding element of each sequence.   

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
This implementation is completely analogous to the previous at the tecnical level, however, instead of iterating on all possible sequences it iterates on a random set. 

#### Categorization by compactness
```julia
function categorize_by_compactness(paths::Vector{Path})::Vector{Vector{Path}}
    self_adjacency_triangles = path_self_adjacency_triangle.(paths)
    compactness = count_adjacent.(self_adjacency_triangles)
    category_count = maximum(compactness)+1
    path_count = length(paths)
    buffer = Matrix{Path}(undef, (category_count, path_count))
    buffer_count = zeros(Int, category_count)
    for i in 1:path_count
        idx = compactness[i]+1
        buffer_count[idx] += 1
        buffer[idx, buffer_count[idx]] = paths[i]
    end
    categories = Vector{Vector{Path}}([buffer[i, 1:buffer_count[i]] for i in 1:category_count])
    return categories
end 
```
This function splits a vector of paths into different compactness "categories", where the compactness of a path is simply calculated as the sum of all elements of the adjacency triangle described above 

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
This function returns the multiplicity frequency of the fundamental state with respect to a set of sequences and one of paths.

### Additional implementations
Here I leave some other implementations that were not practically used in the tests but have been discussed in the theoretical document. In particular both an occupation matrix approach to the extension of self avoiding paths, and, more interestingly, an unoptimized Seed-Sprout approach for the seach of all topological self avoiding paths are presented. More work is needed to assess the efficacy of both methods.

#### State
```julia
mutable struct State
    path::Path
    occupation::Array{Bool}
end
```
A `State` is nothing but the combination of a `Path` and an occupation matrix (`occupation`) this data structure was introduced to implement self avoidance in an efficient manner, however, memory allocation overhead, at least for short sequences turns out to be much more relevant than the cost of the slightly longer naive version. To avoid overloading the memory with potentially useless information we have decided to keep `State` and `Path` as separate structures that are used according to necessities.

#### Seeds search
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

function search_sprouts(length::Int, dim::Int)::Vector{Vector{State}}
    if length > 0
        sprouts_buffer = [Vector{State}(undef, (2*dim)^(i-1)) for i in 1:length]
        sprouts = Vector{Vector{State}}(undef, length)
        sprouts[1] = [State(ZeroPath(dim))]
        for i in 2:length
            count = 0
            for sprout in sprouts[i-1]
                children = extend_SAP(sprout, dim)
                for child in children
                    sprouts_buffer[i][count+1] = child
                    count += 1
                end
            end
            sprouts[i] = sprouts_buffer[i][1:count]
        end
        return sprouts
    else
        return Vector{Vector{State}}([[ZeroState(dim)]])
    end
end

function suture(seed::Path, sprout::Path, dim::Int)
    neck = seed.vertices[seed.length].position
    translated_sprout = path_translate(sprout, neck)
    for seed_vertex in seed.vertices
        for i in 1:Base.length(translated_sprout.vertices)
            if seed_vertex.position == translated_sprout.vertices[i].position && i!=1
                return ZeroPath(dim)
            end
        end
    end
    return path_combine(seed, translated_sprout)
end
```
This is a crude implementation of the seed and sprout searches and of the suture method that is used to connect the two, avoiding self intersection. The peroformance is really promising, however more work would be needed to have a reliable and fast algorithm, so for the moment this is not used in the test scripts. This said it should still work properly, so, if one wants to test it, the funciton `unique_SAPs` contains an implementation of this method that should easily substitute `topological_SAPs`.

### Additional Notes
Considering the scope, efficiency was prioretized over safety, many functions may throw unhandled errors.

At the end I probably just managed to make my inefficient code also unsafe.
