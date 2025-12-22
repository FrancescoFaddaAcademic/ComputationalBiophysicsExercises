
using CairoMakie
using Combinatorics

#Dimension should be smaller than the Sys.WORD_SIZE because of the flag system 
const DIM = 2

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

function ZeroPath()
    return Path{Vertex}([Vertex(fill(0, DIM))],1,UInt(0))
end 

Path(vertices::Vector{Vertex}, asymmetry_flag) = Path{Vertex}(vertices, Base.length(vertices), asymmetry_flag)
Path(vertices::Vector{Vector{Int}}, asymmetry_flag) = Path{Vertex}([Vertex(vertices[i]) for i in 1:length(vertices)], length(vertices), asymmetry_flag)
Path(vertices::Vector{Vector{Int}}, colors::Vector{Int}, asymmetry_flag) = Path{ColoredVertex}([ColoredVertex(vertices[i], colors[i]) for i in 1:length(vertices)], length(vertices), asymmetry_flag)

function State(path)
    occupation = fill(false, fill(path.length*2+1, DIM)...)
    for vertex in path.vertices
        occupation[CartesianIndex((path.length .+ vertex.position .+ 1)...)] = true
    end
    return State(path, occupation)
end

#Check dimension is compatible with the flag system
if DIM >= Sys.WORD_SIZE
    error("The dimension of the system must be smaller than the word size for this processor (=$(Sys.WORD_SIZE))" )
end

function compute_base_vectors(dim::Int)
    base_vectors = Vector{Vector{Int}}(undef, dim)
    for i in 1:dim
        base_vectors[i] = (1:dim) .== i
    end
    return base_vectors
end

function build_orthogonal_matrices()::Vector{Matrix{Int}}
    matrices = Vector{Matrix{Int}}()
    for permutation in permutations(1:DIM)
        for reflection_flag in UInt.(0:2^DIM-1)
            M = zeros(Int, DIM, DIM)
            for dim_idx in 1:DIM
                UInt(2^(dim_idx-1)) & reflection_flag >> dim_idx-1 == UInt(0) ? sign = 1 : sign = -1
                M[dim_idx, permutation[dim_idx]] = sign
            end
            push!(matrices, M)
        end
    end
    return matrices
end

const BASE = compute_base_vectors(DIM) 

const ASYMM_KEY = sum(2 .^(0:DIM-1))

const ORTHO_MATRICES = build_orthogonal_matrices()

#Translates the centered index into that of the gird axis. Bound checking is performed by the Array structure
centered_index(index::Vector{Int}, radius::Int) = radius .+ index .+ 1 

function path_extend(path::Path{Vertex}, step::Vector{Int})::Path
    broken_symmetries_vector = Int.(step .!= 0)
    broken_symmetries_flag = UInt(0)
    for i in eachindex(broken_symmetries_vector)
        broken_symmetries_flag += UInt(broken_symmetries_vector[i] * 2^(i-1))
    end
    return Path(vcat(path.vertices, Vertex(path.vertices[path.length].position .+ step)), path.length + 1, path.asymmetry_flag | broken_symmetries_flag)
end

function path_copy(path::Path)::Path
    return deepcopy(path)
end

function path_translate!(path::Path, translation_vector::Vector{Int})::Path
    for vertex in path.vertices
        vertex.position += translation_vector
    end
    return path
end 

function path_translate(path::Path, translation_vector::Vector{Int})::Path
    return path_translate!(path_copy(path), translation_vector::Vector{Int})
end

function path_concatenate(path1::Path{Vertex}, path2::Path{Vertex})
    last_position = path1.vertices[path1.length].position
    return Path(vcat(path1.vertices, path_translate(path2, last_position).vertices[2:path2.length]), path1.asymmetry_flag | path2.asymmetry_flag)
end

#Extends a self avoiding path 
function extend_SAP(state::State)::Vector{State} 
    path = state.path
    occupation = state.occupation
    new_paths_buffer = Vector{State}(undef, 2*DIM)
    last_position = path.vertices[path.length].position
    buffer_idx = 1
    extensions_count = 0
    for axis_idx in 1:DIM
        for sign in [+1,-1]
            step = sign * BASE[axis_idx]
            new_position = last_position + step
            if !occupation[CartesianIndex(centered_index(new_position, path.length)...)]
                new_paths_buffer[buffer_idx] = State(path_extend(path, step))
                buffer_idx += 1
                extensions_count += 1
            end
        end 
    end
    if extensions_count != 0
        return new_paths_buffer[1:extensions_count]
    else
        return Vector{State}([])
    end
end

function is_completely_asymmetric(path::Path)
    return path.asymmetry_flag == ASYMM_KEY
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

function linear_transform(path::Path, matrix)::Path
    return linear_transform!(deepcopy(path), matrix)
end

function linear_transform!(path::Path, matrix)::Path
    for vertex in path.vertices
        vertex.position = matrix * vertex.position
    end
    return path
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
    for matrix in ORTHO_MATRICES
        if are_equal(long_path, linear_transform(short_path, matrix))
            return true
        end
    end
    return false
end

function filter_equivalent(paths::Vector{Path})::Vector{Path}
    filtered_paths = Vector{Path}()
    for path in paths
        eq_flag = false
        for filtered_path in filtered_paths 
            if are_equivalent(path, filtered_path)
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

function filter_equivalent(states::Vector{State})::Vector{State}
    return State.(filter_equivalent(Vector{Path}([state.path for state in states])))
end

function search_heads(length::Int)::Vector{Vector{State}}
    states = [State(ZeroPath())]            
    heads = [Vector{State}() for _ in 1:length]  
    next_states = Vector{State}()
    next_states_candidates = Vector{State}()

    for i in 1:length-1
        empty!(next_states)
        empty!(next_states_candidates)

        for state in states
            child_states = extend_SAP(state)
            is_first_completely_asymmetric = true
            is_first_symmetry_break = true

            for child_state in child_states
                #MUST BE FIXED
                child_copy = deepcopy(child_state)
                if is_completely_asymmetric(child_copy.path) && is_first_completely_asymmetric
                    push!(heads[i+1], child_copy)
                    is_first_completely_asymmetric = false

                elseif child_copy.path.asymmetry_flag != state.path.asymmetry_flag &&
                       is_first_symmetry_break &&
                       !(is_completely_asymmetric(child_copy.path))

                    push!(next_states, child_copy)
                    is_first_symmetry_break = false

               elseif !((child_state.path.asymmetry_flag != state.path.asymmetry_flag && !is_first_symmetry_break) || (is_completely_asymmetric(child_state.path) && !is_first_completely_asymmetric))
                    push!(next_states_candidates, deepcopy(child_copy))
                end
            end
        end

        append!(next_states, deepcopy.(filter_equivalent(next_states_candidates)))

        states = deepcopy(next_states)

        #println("Heads at iteration $(i+1):")
        #println([head.path for head in heads[i+1]])
        #println("--------------------------------------------------------")
    end
    append!(heads[length], deepcopy(next_states))
    return heads
end

function search_bodies(length::Int)::Vector{Vector{State}}
    bodies = [Vector{State}() for _ in 1:length]
    bodies[1] = [State(ZeroPath())]
    for i in 1:length-1
        for body in bodies[i]
            append!(bodies[i+1], deepcopy(extend_SAP(body)))
        end
    end
    return bodies
end

function unwrap_recursive_vector(hypervector::Vector{Vector{T}})::Vector{T} where (T <: Any) #Bisogna migliorarlo, è un disastro
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

function suture(head::Path, body::Path)
    neck = head.vertices[head.length].position
    translated_body = path_translate(body, neck)
    for head_vertex in head.vertices
        for body_vertex in translated_body.vertices[2:body.length]
            if head_vertex.position == body_vertex.position
                return ZeroPath()
            end
        end
    end
    return path_concatenate(head, body)
end

function suture(head::State, body::State)
    return State(suture(head.path, body.path))
end 

function unchecked_suture(heads::Vector{Vector{State}}, bodies::Vector{Vector{State}}, length::Int)::Vector{Path}
    buffer_size = 0
    for i in 1:length-DIM        
        buffer_size += Base.length(bodies[i]) * Base.length(heads[length-i+1])
    end
    buffer = Vector{Path}(undef, buffer_size)
    buffer_idx = 0
    for i in 1:length-DIM
        for body in bodies[i]
            for head in heads[length-i+1]
                path = deepcopy(suture(head.path, body.path))
                if path.length != 1
                    buffer_idx += 1
                    buffer[buffer_idx] = deepcopy(path)
                end
            end
        end
    end
    return buffer[1:buffer_idx-1]
end

length = 4
heads = search_heads(length)
bodies = search_bodies(length-DIM)
paths = unchecked_suture(heads, bodies, length)