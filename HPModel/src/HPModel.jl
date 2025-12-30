module HPModel

#using BenchmarkTools
#using Profile
using Combinatorics
using StaticArrays

export Vertex, Path, State, search_seeds, search_sprouts, suture, all_combinations_suture, unique_SAPs, all_binary_sequences, random_binary_sequences,
path_self_adjacency_triangle, count_adjacent, calculate_self_interaction, multiplicative_binary_interaction_model, categorize_by_compactness, calculate_g, topological_SAPs

abstract type AbstractVertex end
abstract type AbstractPath end

struct Vertex <: AbstractVertex
    position::Vector{Int}
end

mutable struct ColoredVertex <:AbstractVertex
    position::Vector{Int}
    color::Int
end

mutable struct Path{V<:AbstractVertex}<:AbstractPath
    vertices::Vector{V}
    length::Int
    asymmetry_flag::UInt
end

mutable struct State
    path::Path
    occupation::Array{Bool}
end

#Initializes a trivial path
function ZeroPath(dim::Int)::Path
    return Path{Vertex}([Vertex(fill(0, dim))],1,UInt(0))
end 

#Initializes a trivial state
function ZeroState(dim::Int)::State
    occupation = fill(false, fill(3, dim)...)
    occupation[ntuple(_->2, dim)...] = true
    return State(ZeroPath(dim), occupation)
end

@inline Path(vertices::Vector{Vertex}, asymmetry_flag) = Path{Vertex}(vertices, Base.length(vertices), asymmetry_flag)
@inline Path(vertices::Vector{Vector{Int}}, asymmetry_flag) = Path{Vertex}([Vertex(vertices[i]) for i in 1:Base.length(vertices)], Base.length(vertices), asymmetry_flag)

#Extracts the dimension from the path
function path_dim(path::Path)::Int
    return Base.length(path.vertices[1].position)
end

#Initializes a state from a given path
@inline function State(path::Path)::State
    dim = Base.length(path.vertices[1].position)
    occupation = fill(false, fill(path.length*2+1, dim)...)
    for vertex in path.vertices
        occupation[CartesianIndex((path.length .+ vertex.position .+ 1)...)] = true
    end
    return State(path, occupation)
end

#GEOMETRIC PROPERTIES --------------------------------------------------------------

#Base vectors in a generic dimension
function compute_base_vectors(dim::Int)
    base_vectors = Vector{Vector{Int}}(undef, dim)
    for i in 1:dim
        base_vectors[i] = (1:dim) .== i
    end
    return base_vectors
end

#Orthogonal Matrices in a generic dimension
function O(dim::Int)
    matrices = Vector{Matrix{Int}}()
    for permutation in permutations(1:dim)
        for reflection_flag in UInt.(0:2^dim-1)
            M = zeros(Int, dim, dim)
            for dim_idx in 1:dim
                UInt(2^(dim_idx-1)) & reflection_flag >> dim_idx-1 == UInt(0) ? sign = 1 : sign = -1
                M[dim_idx, permutation[dim_idx]] = sign
            end
            push!(matrices, M)
        end
    end
    mn = matrix_number(dim)
    return SVector{mn, SMatrix{dim, dim, Int}}([SMatrix{dim, dim, Int}(matrices[i]) for i in 1:mn])
end

#Order of the Hyperoctaedral group
matrix_number(dim::Int) = factorial(dim)*2^dim

#Asymmetry key corresponding to full asymmetry
asymm_key(dim::Int) = sum(2 .^(1:dim))

#Translates the centered index into that of the gird axis. Bound checking is performed by the Array structure
@inline centered_index(index::Vector{Int}, radius::Int) = ntuple(i -> radius + index[i] + 1, Base.length(index))

#Extends a path regardless of collisions
function path_extend(path::Path{Vertex}, new_position::Vector{Int}, dim::Integer)::Path
    broken_symmetries_vector = Int.(new_position .!= 0)
    broken_symmetries_flag = UInt(0)
    uniqueness_bit = (path.asymmetry_flag ⊻ asymm_key(dim)) == 0
    for i in eachindex(broken_symmetries_vector)
        broken_symmetries_flag += UInt(broken_symmetries_vector[i] * 2^(i))
    end
    return Path(vcat(path.vertices, Vertex(new_position)), path.length + 1, path.asymmetry_flag | broken_symmetries_flag | uniqueness_bit)
end

#Extends a state regardless of collisions 
function state_extend(state::State, new_position::Vector{Int}, dim::Int)::State
    path = state.path
    len = path.length
    new_len = len + 1 
    cent_index = centered_index(new_position, new_len)
    new_side_length = 2new_len + 1
    new_size = ntuple(_-> new_side_length, dim)
    new_occupation = fill(false, new_size...)
    window = ntuple(_->2:new_side_length - 1, dim)
    new_occupation[window...] = state.occupation
    new_occupation[cent_index...] = true
    new_path = path_extend(path, new_position, dim)
    return State(new_path, new_occupation)
end

#Creates a copy of a path
function path_copy(path::Path)::Path
    return Path(Vector{Vector{Int}}([path.vertices[i].position for i in 1:path.length]), path.asymmetry_flag)
end

#Creates a copy of a state
function state_copy(state::State)::State
    return State(path_copy(state.path), copy(state.occupation))
end

#Translates a path (modifies mutable content!)
function path_translate!(path::Path, translation_vector::Vector{Int})::Path
    for vertex in path.vertices
        vertex.position += translation_vector
    end
    return path
end 

#Creates a translated copy of a path
function path_translate(path::Path, translation_vector::Vector{Int})::Path
    dim = Base.length(translation_vector)
    new_vertices = Vector{Vertex}([Vertex(path.vertices[i].position + translation_vector) for i in 1:path.length])
    return Path(new_vertices, path.asymmetry_flag)
end 

#Creates a path containing all the element of the first argument and all but the first elements of the second
@inline function path_combine(path1::Path, path2::Path)::Path
    return Path(vcat(path1.vertices, path2.vertices[2:path2.length]), path1.asymmetry_flag | path2.asymmetry_flag)
end

#Concatenates tail to seed two paths
function path_concatenate(path1::Path{Vertex}, path2::Path{Vertex})::Path
    last_position = path1.vertices[path1.length].position
    return Path(vcat(path1.vertices, path_translate(path2, last_position).vertices[2:path2.length]), path1.asymmetry_flag | path2.asymmetry_flag)
end

#Extends a self avoiding path using states
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

function is_in(path::Path, position::Vector{Int})::Bool
    for vertex in path.vertices
        if vertex.position == position
            return true
        end
    end
    return false 
end

#Extensd a self avoiding path without states
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

#Asserts if a path has no rigid symmetries
function is_completely_asymmetric(path::Path, dim::Int)
    return path.asymmetry_flag == asymm_key(dim)
end

#Asserts if two paths are equivalent
function are_equal(path1::Path, path2::Path)::Bool
    if path1.length != path2.length
        return false
    end
    length = path1.length
    for i in 1:length
        if path1.vertices[i].position != path2.vertices[i].position
            return false
        end
    end
    return true
end

#Transforms a path accordingly to a specified transformation matrix (modifies mutable content!)
function linear_transform!(path::Path, matrix)::Path
    for vertex in path.vertices
        vertex.position .= matrix * vertex.position
    end
    return path
end

#Creates a copy of a path and transforms it according to a specified transformation matrix
function linear_transform(path::Path, matrix)::Path
    return linear_transform!(deepcopy(path), matrix)
end

#Assertrs if two paths are equivalent up to a set of linear transformations represented by specified matrices 
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

#Asserts if two paths are equivalent up to orthogonal transformations
function are_equivalent(path1::Path, path2::Path)::Bool
    return are_equivalent(path1, path2, path_dim(path1))
end

#Asserts if two paths are equivalent up to orthogonal transformations of specified dimension
function are_equivalent(path1::Path, path2::Path, dim::Int)::Bool
    return are_equivalent(path1, path2, O(dim))::Bool
end

#Fuses avector of vector in a single vector
function fast_flatten(nested::Vector{Vector{T}})::Vector{T} where T 
    len = sum(length, nested)
    flattened = Vector{T}(undef, len)
    offset = 1
    for vector in nested
        flattened[offset:offset+length(vector)-1] = vector
        offset += length(vector)
    end
    return flattened
end

#Selects one single representative for each equivalence class
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

#Constructs all the possible topological self avoiding paths up to a specified length 
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

#ALTERNATIVE METHOD (Functioning but unused) ----------------------------------------------- 

#Constructrs a maximal subset of non-equivalent paths (under orthogonal transformation equivalence) from a given set
function filter_equivalent(paths::Vector{Path}, dim::Int)::Vector{Path}
    filtered_paths = Vector{Path}()
    matrices = O(dim)
    for path in paths
        eq_flag = false
        for filtered_path in filtered_paths 
            if are_equivalent(path, filtered_path, matrices)
                eq_flag = true
                break   
            end
        end
        if eq_flag == false
            push!(filtered_paths, deepcopy(path))
        end
    end
    return filtered_paths
end

#Constructrs a maximal subset of non-equivalent paths (under linear transformations specified by a set of matrices) from a given set
function filter_equivalent(paths::Vector{Path}, dim::Int, matrices)::Vector{Path}
    filtered_paths = Vector{Path}()
    for path in paths
        eq_flag = false
        for filtered_path in filtered_paths 
            if are_equivalent(path, filtered_path, matrices)
                eq_flag = true
                break   
            end
        end
        if eq_flag == false
            push!(filtered_paths, deepcopy(path))
        end
    end
    return filtered_paths
end

#Constructrs a maximal subset of non-equivalent states (under orthogonal transformation equivalence) from a given set
function filter_equivalent(states::Vector{State}, dim::Int)::Vector{State}
    return State.(filter_equivalent(Vector{Path}([state.path for state in states]), dim))
end

#Constructrs a maximal subset of non-equivalent states (under linear transformations specified by a set of matrices) from a given set
function filter_equivalent(states::Vector{State}, dim::Int, matrices)::Vector{State}
    return State.(filter_equivalent(Vector{Path}([state.path for state in states]), dim, matrices))
end

#Finds all the minimal asymmetric paths up to a length together with the symmetric paths of the same length   
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

#Constructs all the self avoiding paths up to a specified length
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

#Combines two paths end to beginning, if there are collisions returns a trivial path
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

#Combines two states end to beginning, if there are collisions returns a trivial state
function suture(seed::State, sprout::State, dim::Int)
    return State(suture(seed.path, sprout.path, dim))
end 

#Combines the elements of two sets of self avoiding paths in all possible ways
function all_combinations_suture(seeds::Vector{Vector{State}}, sprouts::Vector{Vector{State}}, length::Int, dim::Int)::Vector{Path}
    if length > dim
        buffer_size = 0
        for i in 1:length-dim        
            buffer_size += Base.length(sprouts[i]) * Base.length(seeds[length-i+1])
        end
        buffer = Vector{Path}(undef, buffer_size)
        buffer_idx = 0
        for i in 1:length-dim
            for sprout in sprouts[i]
                for seed in seeds[length-i+1]
                    path = suture(seed.path, sprout.path, dim)
                    if path.length != 1
                        buffer_idx += 1
                        buffer[buffer_idx] = path
                    end
                end
            end
        end
        return buffer[1:buffer_idx]
    else 
        return [seeds[dim][i].path for i in 1:Base.length(seeds[dim])]
    end
end

#Construct all the possible self avoiding paths (modulo orthogonal group) up to a length
function unique_SAPs(depth::Int, dim::Int)
    seeds = search_seeds(depth, dim)
    sprouts = search_sprouts(depth-dim, dim)
    all_combinations = all_combinations_suture(seeds, sprouts, depth, dim)
    return all_combinations
end

#SELF INTERACTION METHODS ----------------------------------------------------------------------

#Evaluates norm one distance between vertices
function one_distance(vertex1::Vertex, vertex2::Vertex, dim::Integer)::Integer
    distance = 0
    for i in 1:dim
        distance += abs(vertex1.position[i] - vertex2.position[i])
    end
    return distance
end

#Constructs of the self adjacency triangle
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

#Calculates the number of self adjacencies of a path by providing its adjacency matrix 
function count_adjacent(triangle::Vector{Vector{Bool}})
    count = 0
    for i in 1:Base.length(triangle)
        for j in 1:i
            count += triangle[i][j]
        end
    end
    return count
end

#Evaluates the cumulative effect of a first neighbor self interaction with respect to a given interaction model 
function calculate_self_interaction(sequence::Vector{Int}, triangle::Vector{Vector{Bool}}, interaction_model)
    interaction = 0
    for i in 1:Base.length(triangle)
        for j in 1:i
            interaction += triangle[i][j] * interaction_model(sequence[i], sequence[j])
        end
    end
    return interaction
end

#INTERACTION MODELS --------------------------------------------------------------------

#Simple multiplicative interaction
function multiplicative_binary_interaction_model(a,b)
    return a*b
end

#Interaction model by lee er al.
function lee_binary_interaction(a,b)
    if a == b == 1
        return 2.3
    elseif a == b == 0
        return 0
    else
        return 1
    end
end

#SEQUENCES METHODS ------------------------------------------------------------

#Constructs all binary sequences of a specified lenght
function all_binary_sequences(len::Int)::Vector{Vector{Int}}
    sequences = Vector{Vector{Int}}([zeros(len) for _ in 1:2^len])
    for i in 1:2^len
        for j in 1:len
            sequences[i][j] = Int(UInt(i) >> (j-1) & UInt(1))
        end
    end
    return sequences
end

#Constructs a random vector of specified length whose elements are binary sequences of a specified length
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

#Organizes paths by compactness
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

#Calculates the number of sequences corresponding to any multiplicity of configurations of maximal interaction
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

#Calculates the number of sequences corresponding to any multiplicity of configurations of maximal interaction
function calculate_g(sequences::Vector{Vector{Int}}, paths::Vector{Path}, interaction_model)
    return calculate_g(sequences, path_self_adjacency_triangle.(paths), interaction_model)
end

#UTILS -----------------------------------------------------------------------------

#Specialized printing function for path visualization
function print_path(path::Path)
    print("<")
    for vertex in path.vertices
        print(vertex.position)
    end
    print(">\n")
end

#Specialized printing function for adjacency visualization
function print_triangle(triangle::Vector{Vector{Bool}})
    for i in 1:length(triangle)
        for j in 1:i
            triangle[i][j] ? print("X") : print("O")
        end
        print("\n")
    end
end

end