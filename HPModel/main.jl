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

function ZeroPath(dim::Int)::Path
    return Path{Vertex}([Vertex(fill(0, dim))],1,UInt(0))
end 

function ZeroState(dim::Int)::State
    occupation = fill(false, fill(3, dim)...)
    occupation[ntuple(_->2, dim)...] = true
    return State(ZeroPath(dim), occupation)
end

Path(vertices::Vector{Vertex}, asymmetry_flag) = Path{Vertex}(vertices, Base.length(vertices), asymmetry_flag)
Path(vertices::Vector{Vector{Int}}, asymmetry_flag) = Path{Vertex}([Vertex(vertices[i]) for i in 1:Base.length(vertices)], Base.length(vertices), asymmetry_flag)
Path(vertices::Vector{Vector{Int}}, colors::Vector{Int}, asymmetry_flag) = Path{ColoredVertex}([ColoredVertex(vertices[i], colors[i]) for i in 1:length(vertices)], length(vertices), asymmetry_flag)

function path_dim(path::Path)::Int
    return Base.length(path.vertices[1].position)
end

function State(path::Path)::State
    dim = Base.length(path.vertices[1].position)
    occupation = fill(false, fill(path.length*2+1, dim)...)
    for vertex in path.vertices
        occupation[CartesianIndex((path.length .+ vertex.position .+ 1)...)] = true
    end
    return State(path, occupation)
end

#Geometric properties--------------------------------------------------------------

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


asymm_key(dim::Int) = sum(2 .^(0:dim-1))

#Translates the centered index into that of the gird axis. Bound checking is performed by the Array structure
centered_index(index::Vector{Int}, radius::Int) = ntuple(i -> radius + index[i] + 1, Base.length(index))

function path_extend(path::Path{Vertex}, new_position::Vector{Int})::Path
    broken_symmetries_vector = Int.(new_position .!= 0)
    broken_symmetries_flag = UInt(0)
    for i in eachindex(broken_symmetries_vector)
        broken_symmetries_flag += UInt(broken_symmetries_vector[i] * 2^(i-1))
    end
    return Path(vcat(path.vertices, Vertex(new_position)), path.length + 1, path.asymmetry_flag | broken_symmetries_flag)
end

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

function path_copy(path::Path)::Path
    return Path(Vector{Vector{Int}}([path.vertices[i].position for i in 1:path.length]), path.asymmetry_flag)
end

function state_copy(state::State)::State
    return State(path_copy(state.path), copy(state.occupation))
end

function path_translate!(path::Path, translation_vector::Vector{Int})::Path
    for vertex in path.vertices
        vertex.position += translation_vector
    end
    return path
end 

function path_translate(path::Path, translation_vector::Vector{Int})::Path
    dim = Base.length(translation_vector)
    new_vertices = Vector{Vertex}([Vertex(path.vertices[i].position + translation_vector) for i in 1:path.length])
    return Path(new_vertices, path.asymmetry_flag)
end 

function path_add(path1::Path, path2::Path)::Path
    return Path(vcat(path1.vertices, path2.vertices[2:path2.length]), path1.asymmetry_flag | path2.asymmetry_flag)
end

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
    if extensions_count != 0
        return new_state_buffer[1:extensions_count]
    else
        return Vector{State}([])
    end
end

function is_completely_asymmetric(path::Path, dim::Int)
    return path.asymmetry_flag == asymm_key(dim)
end

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

function linear_transform!(path::Path, matrix)::Path
    for vertex in path.vertices
        vertex.position = matrix * vertex.position
    end
    return path
end

function linear_transform(path::Path, matrix)::Path
    return linear_transform!(path_copy(path), matrix)
end

function are_equivalent(path1::Path, path2::Path)::Bool
    #Order paths by length for higher efficiency
    if path1.length < path2.length
        short_path = path1
        long_path = path2
    else
        short_path = path2
        long_path = path1
    end
    for matrix in O(dim(path_dim(path1)))
        if are_equal(long_path, linear_transform(short_path, matrix))
            return true
        end
    end
    return false
end

function are_equivalent(path1::Path, path2::Path, dim::Int)::Bool
    #Order paths by length for higher efficiency
    if path1.length < path2.length
        short_path = path1
        long_path = path2
    else
        short_path = path2
        long_path = path1
    end
    for matrix in O(dim)
        if are_equal(long_path, linear_transform(short_path, matrix))
            return true
        end
    end
    return false
end

function are_equivalent(path1::Path, path2::Path, matrices)::Bool
    #Order paths by length for higher efficiency
    if path1.length < path2.length
        short_path = path1
        long_path = path2
    else
        short_path = path2
        long_path = path1
    end
    for matrix in matrices
        if are_equal(long_path, linear_transform(short_path, matrix))
            return true
        end
    end
    return false
end

#MUST BE IMPROVED
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

function filter_equivalent(states::Vector{State}, dim::Int)::Vector{State}
    return State.(filter_equivalent(Vector{Path}([state.path for state in states]), dim))
end

function filter_equivalent(states::Vector{State}, dim::Int, matrices)::Vector{State}
    return State.(filter_equivalent(Vector{Path}([state.path for state in states]), dim, matrices))
end

function search_heads(length::Int, dim::Int)::Vector{Vector{State}}          
    heads = Vector{Vector{State}}(undef, length) 
    if dim > 2
        max_size = Base.length(unique_SAPs(length,dim-1))*2
    else
        max_size = 1*2
    end
    matrices = O(dim)
    heads_buffer = Vector{State}(undef, max_size)
    states = Vector{State}(undef, max_size)
 
    next_states = Vector{State}(undef, max_size)
    next_states_candidates = Vector{State}(undef, max_size)

    states[1] = ZeroState(dim)
    states_count = 1

    heads_count = 0
    next_states_count = 0
    for i in 1:length-1
        heads[i] = heads_buffer[1:heads_count]
        next_states_count = 0
        next_states_candidates_count = 0
        heads_count = 0
        for j in 1:states_count
            parent = states[j]
            children = extend_SAP(parent, dim)
            is_first_completely_asymmetric = true
            is_first_symmetry_break = true
            for child in children
                ica = is_completely_asymmetric(child.path, dim)
                symmetry_break = child.path.asymmetry_flag != parent.path.asymmetry_flag
                if ica && is_first_completely_asymmetric
                    heads_count += 1
                    heads_buffer[heads_count] = state_copy(child)
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
    heads_buffer[heads_count+1:heads_count+next_states_count] = next_states[1:next_states_count]
    heads_count += next_states_count
    heads[length] = heads_buffer[1:heads_count]
    return heads
end


function search_bodies(length::Int, dim::Int)::Vector{Vector{State}}
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
end

#DA SISTEMARE
function unwrap_recursive_vector(hypervector::Vector{Vector{T}})::Vector{T} where (T <: Any) 
    out_len = 0
    for vector in hypervector
       out_len += size(vector, 1)
    end
    out_vector = Vector{T}(undef, out_len)
    offset = 0
    for vector in hypervector
        for i in 1:size(vector, 1)
            out_vector[offset + i] = vector[i]
        end
        offset += size(vector, 1)
    end
    return out_vector
end

function suture(head::Path, body::Path, dim::Int)
    neck = head.vertices[head.length].position
    translated_body = path_translate(body, neck)
    for head_vertex in head.vertices
        for i in 1:Base.length(translated_body.vertices)
            if head_vertex.position == translated_body.vertices[i].position && i!=1
                return ZeroPath(dim)
            end
        end
    end
    return path_add(head, translated_body)
end

function suture(head::State, body::State, dim::Int)
    return State(suture(head.path, body.path, dim))
end 

function unchecked_suture(heads::Vector{Vector{State}}, bodies::Vector{Vector{State}}, length::Int, dim::Int)::Vector{Path}
    buffer_size = 0
    for i in 1:length-dim        
        buffer_size += Base.length(bodies[i]) * Base.length(heads[length-i+1])
    end
    buffer = Vector{Path}(undef, buffer_size)
    buffer_idx = 0
    for i in 1:length-dim
        for body in bodies[i]
            for head in heads[length-i+1]
                path = suture(head.path, body.path, dim)
                if path.length != 1
                    buffer_idx += 1
                    buffer[buffer_idx] = path
                end
            end
        end
    end
    return buffer[1:buffer_idx]
end

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

function print_triangle(triangle::Vector{Vector{Bool}})
    for i in 1:length(triangle)
        for j in 1:i
            if triangle[i][j]
                print("X")
            else
                print("O")
            end
        end
        print("\n")
    end
end

function count_adjacent(triangle::Vector{Vector{Bool}})
    count = 0
    for i in 1:Base.length(triangle)
        for j in 1:i
            if triangle[i][j]
                count += 1
            end
        end
    end
    return count
end

function calculate_self_interaction(sequence::Vector{Int}, triangle::Vector{Vector{Bool}}, interaction_model)
    interaction = 0
    for i in 1:Base.length(triangle)
        for j in 1:i
            if triangle[i][j]
                interaction += interaction_model(sequence[i], sequence[j])
            end
        end
    end
    return interaction
end

function print_path(path::Path)
    print("<")
    for vertex in path.vertices
        print(vertex.position)
    end
    print(">\n")
end

function unique_SAPs(depth::Int, dim::Int)
    heads = search_heads(depth, dim)
    bodies = search_bodies(depth-dim, dim)
    return unchecked_suture(heads, bodies, depth, dim)
end

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

function categorize_by_compactedness(SAPs::Vector{Path})::Vector{Vector{Path}}
    self_adjacency_triangles = path_self_adjacency_triangle.(SAPs)
    adjacency_count = count_adjacent.(self_adjacency_triangles)
    categories = Vector{Vector{Path}}([[] for _ in 1:maximum(adjacency_count)+1])
    for i in 1:Base.length(SAPs)
        push!(categories[adjacency_count[i]+1], path_copy(SAPs[i]))
    end
    return categories
end 

function plot_compactedness_categorized(categorized_paths::Vector{Vector{Path}})
    bins = Base.length(categorized_paths)
    fig = GLMakie.Figure(size=(800, 600), fontsize = 25)
    ax = Axis(
        fig[1, 1], 
        title = "Number of sequences by compactedness", 
        xlabel = L"\rho", 
        ylabel = "Number of sequences",
        xgridvisible=false,
        ygridvisible=false,
        xticklabelsize = 20,
        yticklabelsize = 15,
        xticks = round.(0:1/(bins-1):1, digits = 2),
        yticklabelrotation = π/4,
        )
    xlims!(ax, -0.1, 1.1)
    ylims!(ax, 0, nothing)
    barplot!(ax, (0:bins-1)./(bins-1), Base.length.(categorized_paths))
    return fig
end

function plot_compactedness(paths::Vector{Path})
    categorized_paths = categorize_by_compactedness(paths)
    plot_compactedness_categorized(categorized_paths)
end

len = 16
dim = 2

SAPs = unique_SAPs(len, dim)
categorized_SAPs = categorize_by_compactedness(SAPs)

plot_compactedness_categorized(categorized_SAPs)

self_adjacency_triangles = path_self_adjacency_triangle.(categorized_SAPs[Base.length(categorized_SAPs)])

sequences = all_binary_sequences(len)

function interaction_model(a,b)
    return a*b
end

g = zeros(Int, Base.length(categorized_SAPs[Base.length(categorized_SAPs)]))

count_unique_native = 0
for sequence in sequences
    Emax = 0
    multiplicity = 0
    for i in 1:Base.length(categorized_SAPs[Base.length(categorized_SAPs)])
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

fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1])
xlims!(ax, 0, 21)
barplot!(ax, g)
