module HPModel

using BenchmarkTools
using Profile
using GLMakie
using CairoMakie
using Combinatorics
using StaticArrays
using LaTeXStrings

abstract type AbstractVertex end
abstract type AbstractPath end

mutable struct Vertex <: AbstractVertex
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

Path(vertices::Vector{Vertex}, asymmetry_flag) = Path{Vertex}(vertices, Base.length(vertices), asymmetry_flag)
Path(vertices::Vector{Vector{Int}}, asymmetry_flag) = Path{Vertex}([Vertex(vertices[i]) for i in 1:Base.length(vertices)], Base.length(vertices), asymmetry_flag)
Path(vertices::Vector{Vector{Int}}, colors::Vector{Int}, asymmetry_flag) = Path{ColoredVertex}([ColoredVertex(vertices[i], colors[i]) for i in 1:length(vertices)], length(vertices), asymmetry_flag)

#Extracts the dimension from the path
function path_dim(path::Path)::Int
    return Base.length(path.vertices[1].position)
end

#Initializes a state from a given path
function State(path::Path)::State
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
asymm_key(dim::Int) = sum(2 .^(0:dim-1))

#Translates the centered index into that of the gird axis. Bound checking is performed by the Array structure
centered_index(index::Vector{Int}, radius::Int) = ntuple(i -> radius + index[i] + 1, Base.length(index))

#Extends a path regardless of collisions
function path_extend(path::Path{Vertex}, new_position::Vector{Int})::Path
    broken_symmetries_vector = Int.(new_position .!= 0)
    broken_symmetries_flag = UInt(0)
    for i in eachindex(broken_symmetries_vector)
        broken_symmetries_flag += UInt(broken_symmetries_vector[i] * 2^(i-1))
    end
    return Path(vcat(path.vertices, Vertex(new_position)), path.length + 1, path.asymmetry_flag | broken_symmetries_flag)
end

#Extends a state regardless of collisions 
function state_extend(state::State, new_position::Vector{Int}, dim::Int)::State
    path = state.path
    new_side_length = (path.length+1)*2+1
    new_size = ntuple(_-> new_side_length, dim)
    new_occupation = fill(false, new_size...)
    window = ntuple(_->2:(path.length*2+1)+1, dim)
    new_occupation[window...] = state.occupation
    new_occupation[centered_index(new_position, path.length+1)...] = true
    new_path = path_extend(path, new_position)
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
function path_combine(path1::Path, path2::Path)::Path
    return Path(vcat(path1.vertices, path2.vertices[2:path2.length]), path1.asymmetry_flag | path2.asymmetry_flag)
end

#Concatenates tail to seed two paths
function path_concatenate(path1::Path{Vertex}, path2::Path{Vertex})::Path
    last_position = path1.vertices[path1.length].position
    return Path(vcat(path1.vertices, path_translate(path2, last_position).vertices[2:path2.length]), path1.asymmetry_flag | path2.asymmetry_flag)
end

#Extends a self avoiding path 
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
        vertex.position = matrix * vertex.position
    end
    return path
end

#Creates a copy of a path and transforms it according to a specified transformation matrix
function linear_transform(path::Path, matrix)::Path
    return linear_transform!(path_copy(path), matrix)
end

#Assertrs if two paths are equivalent up to a set of linear transformations represented by specified matrices 
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

#Asserts if two paths are equivalent up to orthogonal transformations
function are_equivalent(path1::Path, path2::Path)::Bool
    return are_equivalent(path1, path2, path_dim(path1))
end

#Asserts if two paths are equivalent up to orthogonal transformations of specified dimension
function are_equivalent(path1::Path, path2::Path, dim::Int)::Bool
    return are_equivalent(path1, path2, O(dim))::Bool
end



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
function search_bodies(length::Int, dim::Int)::Vector{Vector{State}}
    if length > 0
        bodies_buffer = [Vector{State}(undef, (2*dim)^(i-1)) for i in 1:length]
        bodies = Vector{Vector{State}}(undef, length)
        bodies[1] = [State(ZeroPath(dim))]
        for i in 2:length
            count = 0
            for body in bodies[i-1]
                children = extend_SAP(body, dim)
                for child in children
                    bodies_buffer[i][count+1] = child
                    count += 1
                end
            end
            bodies[i] = bodies_buffer[i][1:count]
        end
        return bodies
    else
        return Vector{Vector{State}}([[ZeroState(dim)]])
    end
end

#Combines two paths end to beginning, if there are collisions returns a trivial path
function suture(seed::Path, body::Path, dim::Int)
    neck = seed.vertices[seed.length].position
    translated_body = path_translate(body, neck)
    for seed_vertex in seed.vertices
        for i in 1:Base.length(translated_body.vertices)
            if seed_vertex.position == translated_body.vertices[i].position && i!=1
                return ZeroPath(dim)
            end
        end
    end
    return path_combine(seed, translated_body)
end

#Combines two states end to beginning, if there are collisions returns a trivial state
function suture(seed::State, body::State, dim::Int)
    return State(suture(seed.path, body.path, dim))
end 

#Combines the elements of two sets of self avoiding paths in all possible ways
function all_combination_suture(seeds::Vector{Vector{State}}, bodies::Vector{Vector{State}}, length::Int, dim::Int)::Vector{Path}
    if length > dim
        buffer_size = 0
        for i in 1:length-dim        
            buffer_size += Base.length(bodies[i]) * Base.length(seeds[length-i+1])
        end
        buffer = Vector{Path}(undef, buffer_size)
        buffer_idx = 0
        for i in 1:length-dim
            for body in bodies[i]
                for seed in seeds[length-i+1]
                    path = suture(seed.path, body.path, dim)
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
    bodies = search_bodies(depth-dim, dim)
    all_combinations = all_combination_suture(seeds, bodies, depth, dim)
    return all_combinations
end

#SELF INTERACTION METHODS ----------------------------------------------------------------------

#Calculates the upper triangle of the adjacency matrix of a path 
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

#SEQUENCES METHODS ------------------------------------------------------------

#Constructs all binary sequences of a specified lenght
function all_binary_sequences(len::Int)::Vector{Vector{Int}}
    sequences = Vector{Vector{Int}}([zeros(len) for _ in 1:2^len])
    print(typeof(sequences))
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

#Organizes paths by compactedness
function categorize_by_compactedness(paths::Vector{Path})::Vector{Vector{Path}}
    self_adjacency_triangles = path_self_adjacency_triangle.(paths)
    adjacency_count = count_adjacent.(self_adjacency_triangles)
    categories = Vector{Vector{Path}}([[] for _ in 1:maximum(adjacency_count)+1])
    for i in 1:Base.length(paths)
        push!(categories[adjacency_count[i]+1], path_copy(paths[i]))
    end
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


function plot_compactedness_categorized(categorized_paths::Vector{Vector{Path}})
    bins = Base.length(categorized_paths)
    binwidth = 1/(bins-1)
    fig = GLMakie.Figure(size=(800, 600), fontsize = 25)
    ax = Axis(
        fig[1, 1], 
        #title = "Number of sequences by compactedness", 
        xlabel = L"\rho", 
        ylabel = "Number of sequences",
        xgridvisible=false,
        ygridvisible=true,
        xticklabelsize = 20,
        yticklabelsize = 18,
        xticks = round.(0:1/(bins-1):1, digits = 2),
        yticklabelrotation = π/4,
        )
    xlims!(ax, -binwidth/2, 1+binwidth/2)
    ylims!(ax, 0, nothing)
    barplot!(ax, (0:bins-1)./(bins-1), Base.length.(categorized_paths))
    return fig
end

function plot_compactedness(paths::Vector{Path})
    categorized_paths = categorize_by_compactedness(paths)
    plot_compactedness_categorized(categorized_paths)
end

function plot_g(g::Vector{Int})
    bins = Base.length(g)
    binwidth = 1
    fig = Figure(size=(800, 600), fontsize = 25)
    ax = Axis(
        fig[1, 1],
        xlabel = L"g(s)", 
        ylabel = "Number of sequences",
        xgridvisible=false,
        ygridvisible=true,
        xticklabelsize = 20,
        yticklabelsize = 20,
        xticks = round.(1:ceil(bins/15):bins, digits = 2),
    )
    xlims!(ax, 1-binwidth/2, bins+binwidth/2)
    ylims!(ax, 0, nothing)
    barplot!(ax, g)
    return fig
end

len = 10
dim = 2

SAPs = unique_SAPs(len, dim)
categorized_SAPs = categorize_by_compactedness(SAPs)

plot_compactedness_categorized(categorized_SAPs)

self_adjacency_triangles = path_self_adjacency_triangle.(SAPs)

sequences = all_binary_sequences(len)
sequences = random_binary_sequences(len, 10000)
g = calculate_g(sequences, self_adjacency_triangles, multiplicative_binary_interaction_model)
plot_g(g)

end