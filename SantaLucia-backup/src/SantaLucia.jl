module SantaLucia

using StaticArrays
using CSV
using DataFrames


const word_size = Sys.WORD_SIZE
const bases = SVector{4,Char}(['a', 'g', 'c', 't'])
const keys = SVector{4,UInt}([0, 1, 2, 3])

#Si noti che bases[keys[i]+1] restituisce bases[i]
const masks = SVector{32,UInt}([2^(2*i)+2^(2*i + 1) for i in 0:(word_size/2)-1])

mutable struct Segment
    words::Vector{UInt}
    length::UInt
    granularity::UInt
end

struct Sequence
    binary::Segment
end

struct PairSequence 
    literal::Vector{String}
    binary::Segment
    length::UInt
end

struct SingleStrand
    sequence::Vector{UInt}
    structure::Vector{UInt}
    length::UInt
end

struct DoubleStrand
    pair_sequence::Vector{UInt}
    length::UInt
end

function ZeroSegment(length, granularity)::Segment
    return Segment(zeros(UInt, ceil(UInt, length * granularity / word_size)), length, granularity)
end

ZeroSegment1(length) = ZeroSegment(length, 2)
ZeroSegment2(length) = ZeroSegment(length, 4)

#Bitwise method to put a key into a word in a certain position
function pack(word::UInt, key::UInt, position::UInt, keysize::Int)::UInt
    return key << (keysize*(position-1)) | (word & ~masks[position])
end

pack(word, key, position, keysize) = pack(UInt(word), UInt(key), UInt(position), Int(keysize))
pack1(word, key, position) = pack(word, key, position, 2)
pack2(word, key, position) = pack(word, key, position, 4)

#Bitwise method to unpack a key out of a word form a certain position
function unpack(word::UInt, position::UInt, keysize::UInt)::UInt
    return (word & masks[position]) >> (keysize*(position-1))
end

#unpack(word, position, keysize) = invokelatest(unpack(UInt(word), UInt(position), Int(keysize)))
unpack1(word, position) = unpack(word, position, 2)
unpack2(word, position) = unpack(word, position, 4)

function segment_outer_index(segment::Segment, position::UInt)
    word_capacity = div(word_size, segment.granularity)
    return cld(position, word_capacity)
end

function segment_inner_index(segment::Segment, position::UInt)
    word_capacity = div(word_size, segment.granularity)
    return mod1(position, word_capacity)
end

function segment_index(segment::Segment, position::UInt)::Tuple{UInt, UInt}
    word_capacity = div(word_size, segment.granularity)
    outer_idx = cld(position, word_capacity)
    inner_idx = mod1(position, word_capacity)
    return (outer_idx, inner_idx)
end

#Estensione dei precedenti metodi a segmenti
function pack!(segment::Segment, key::UInt, position::UInt)::Segment
    (outer_idx, inner_idx) = segment_index(segment, position)
    segment.words[outer_idx] = pack(segment.words[outer_idx], key, inner_idx, segment.granularity)
    return segment
end

pack!(segment, key, position) = pack!(segment, UInt(key), UInt(position))
pack(segment, key, position) = pack!(copy(segment), key, position)

function unpack(segment::Segment, position::UInt)::UInt
    (outer_idx, inner_idx) = segment_index(segment, position)
    return unpack(segment.words[outer_idx], inner_idx, segment.granularity)
end

#Constructs the key corresoinding to the symbol
function key1(symbol::Char)::UInt
    for i in 1:4
        if symbol == bases[i]
            return UInt(i-1)
        end
    end
    error("Unrecognized symbol")
end

#Construct the symbol corresponding to the key
function symbol1(key::UInt)::Char
    return bases[key+1]
end

function Sequence1(string::String)::Sequence
    len = length(string)
    segment::Segment = ZeroSegment(len, 2)
    for i in 1:len
        pack!(segment, key1(string[i]), i)
    end
    return Sequence(segment)
end

function rotate(segment::Segment)::Segment
    rotated_segment = ZeroSegment(segment.length, segment.granularity)
    for i in 1:segment.length
        pack!(rotated_segment, unpack(segment, i), segment.length-i+1)
    end
    return rotated_segment
end

function rotate(sequence::Sequence)::Sequence
    return Sequence(rotate(sequence.binary))
end

function complementary_segment(segment::Segment)::Segment
    Segment([~word for word in segment.words], segment.length, segment.granularity)
end

function complementary_sequence(sequence::Sequence)::Sequence
    return Sequence(complementary_segment(sequence.binary))
end

function bitwise_binary_check(f::Function, segment1::Segment, segment2::Segment)::Bool
    if segment1.length != segment2.length
        return false
    end
    word_capience = div(word_size, segment1.granularity)
    head_size = rem(segment1.length, word_capience)
    word_count = length(segment1.words)
    for i in 1:word_count-1
        if f(segment1.words[i], segment2.words[i]) != 0
            return false
        end
    end
    if f(segment1.words[word_count], segment2.words[word_count]) << (word_size - head_size) != 0
        return false
    end
    return true
end

function bitwise_binary_check(f::Function, sequence1::Sequence, sequence2::Sequence)::Bool
    return bitwise_binary_check(f, sequence1.binary, sequence2.binary)
end

are_equal(s1, s2) = bitwise_binary_check((a,b)->(a ⊻ b), s1, s2)

are_complementary(s1, s2) = bitwise_binary_check((a,b)->~(a ⊻ b), s1, rotate(s2))

function calculate_inner_potentials(sequence::Sequence, data::DataFrame)::Vector{Float64}
    segment = sequence.binary
    ΔH = 0.0
    ΔS = 0.0
    for i in 1:segment.length-1
        term_key = unpack(segment, i) | (unpack(segment, i+1) << segment.granularity)
        ΔH += data[term_key+1,"DeltaH"]
        ΔS += data[term_key+1,"DeltaS"]
    end
    return [ΔH, ΔS]
end

function calculate_terminal_AT_penalty(sequence::Sequence)::Vector{Float64}
    segment = sequence.binary
    len = segment.length
    term_key = unpack(segment, len)
    if term_key == 0
        return [2.2, 6.9]
    end
    return [0.0, 0.0]
end

function calculate_symmetry_correction(sequence::Sequence)::Vector{Float64}
    if are_complementary(sequence, sequence)
        return [0.0, -1.4]
    end
    return [0.0, 0.0]
end

function calculate_thermodynamic_potentials(sequence::Sequence, data_inner::DataFrame)::Vector{Float64}
    Δ = [0.2, -5.7]
    Δ += calculate_inner_potentials(sequence, data_inner)
    Δ += calculate_terminal_AT_penalty(sequence)
    Δ += calculate_symmetry_correction(sequence)
    return Δ
end

function meltingTemperature(H, S, C)
    return H/(S+R*log(C/4))
end

function meltingCurve(β, H, S)::Float64
    x = exp(-(β*H - S)/R)
    return 1-2(-1+sqrt(1+2*x))/(2*x)
end

sequence = Sequence1("cgttga") 

data_filepath = "SL_Inner_Potential.csv"
data_inner = CSV.read(data_filepath, DataFrame)
calculate_inner_potentials(sequence, data_inner)

potentials = calculate_thermodynamic_potentials(sequence, data_inner)

complementary = complementary_sequence(sequence)

df = CSV.read("SL_Parameters.csv", DataFrame)

df = DataFrame(
    term = [String(symbol1.([unpack1(UInt(i), UInt(1)), unpack1(UInt(i), UInt(2))])) for i in 0:15],
    deltaH = [Float64(0) for _ in 0:15],
    deltaS = [Float64(0) for _ in 0:15],
)

CSV.write("SL_Parameters.csv", df)

end