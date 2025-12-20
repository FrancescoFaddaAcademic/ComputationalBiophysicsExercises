
abstract type Vertex end

struct RawVertex
    position::Vector{Int32}
end

struct ColouredVertex
    is_hydrophilic::Bool
    position::Vector{Int32}
end

struct Path{V<:Vertex}
    vertices::Vector{V}
end

const DIM = 3

function compute_base_vectors(dim::Int)
    base_vectors = Vector{Vector{Int}}(undef, dim)
    for i in 1:dim
        base_vectors[i] = (1:dim) .== i
    end
    return base_vectors
end

const BASE_VECS = compute_base_vectors(DIM) 

length = 10

protein = Vector{Vertex}(undef, length)

grid_size::UInt32 = length * 2 + 1 

centered_index(index::Int) = length + index + 1 #Bound checking performed by the Array structure



function raw_paths_up_to_length(length::Integer) #Gets all the (topological) configurations
    raw_paths = Vector{Path{RawVertex}}(undef, length)
    if length == 2
        return [Path{RawVertex}([ZERO_VEC, ZERO_VEC .+ BASE_VECS[1]])]
    end
end

function adjacentcy_matrix(path::Path)

end