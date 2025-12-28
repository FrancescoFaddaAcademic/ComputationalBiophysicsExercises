module SantaLucia

using StaticArrays
using CSV
using DataFrames

export 
Segment,
Sequence,
print_segment,
print_sequence,
invert, 
complementary_segment, 
complementary_sequence,
bitwise_binary_check, 
inner_potentials, 
terminal_AT_penalty, 
symmetry_correction, 
thermodynamic_potentials, 
melting_temperature, 
melting_curve, 
create_NN_template

const word_size = Sys.WORD_SIZE
const bases_symbols = SVector{4,Char}(['A', 'G', 'C', 'T'])
const dotparen_symbols = SVector{3, Char}([')','.','('])
const keys = SVector{4,UInt}([0, 1, 2, 3])

#Si noti che bases[keys[i]+1] restituisce bases[i]
const masks = SVector{32,UInt}([2^(2*i)+2^(2*i + 1) for i in 0:(word_size/2)-1])

const R = 1.987

mutable struct Segment
    words::Vector{UInt}
    length::UInt
    granularity::UInt
end

mutable struct Sequence
    binary::Segment
end

mutable struct Structure
    binary::Segment
end

mutable struct HairPin
    loop_length::UInt
    paired_sequence::Sequence
end

mutable struct SingleStrand
    sequence::Sequence
    structure::Structure
    length::UInt
end

mutable struct DoubleStrand
    sequence::Segment
    length::UInt
end

function ZeroSegment(length, granularity)::Segment
    return Segment(zeros(UInt, ceil(UInt, length * granularity / word_size)), length, granularity)
end

ZeroSegment2(length) = ZeroSegment(length, 2)
ZeroSegment4(length) = ZeroSegment(length, 4)

#Bitwise method to put a key into a word in a certain position
function pack(word::UInt, key::UInt, position::Integer, keysize::Integer)::UInt
    return key << (keysize*(position-1)) | (word & ~masks[position])
end

pack(word, key, position, keysize) = pack(UInt(word), UInt(key), UInt(position), Int(keysize))
pack2(word, key, position) = pack(word, key, position, 2)
pack4(word, key, position) = pack(word, key, position, 4)

#Bitwise method to unpack a key out of a word form a certain position
function unpack(word::UInt, position::Integer, keysize::Integer)::UInt
    return (word & masks[position]) >> (keysize*(position-1))
end

#unpack(word, position, keysize) = invokelatest(unpack(UInt(word), UInt(position), Int(keysize)))
unpack2(word, position) = unpack(word, position, 2)
unpack4(word, position) = unpack(word, position, 4)

function segment_outer_index(segment::Segment, position::Integer)
    word_capacity = div(word_size, segment.granularity)
    return cld(position, word_capacity)
end

function segment_inner_index(segment::Segment, position::Integer)
    word_capacity = div(word_size, segment.granularity)
    return mod1(position, word_capacity)
end

function segment_index(segment::Segment, position::Integer)::Tuple{UInt, UInt}
    word_capacity = div(word_size, segment.granularity)
    outer_idx = cld(position, word_capacity)
    inner_idx = mod1(position, word_capacity)
    return (outer_idx, inner_idx)
end

#Estensione dei precedenti metodi a segmenti
function pack!(segment::Segment, key::UInt, position::Integer)::Segment
    (outer_idx, inner_idx) = segment_index(segment, position)
    segment.words[outer_idx] = pack(segment.words[outer_idx], key, inner_idx, segment.granularity)
    return segment
end

pack!(segment, key, position) = pack!(segment, UInt(key), UInt(position))
pack(segment, key, position) = pack!(deepcopy(segment), key, position)

function unpack(segment::Segment, position::Integer)::UInt
    (outer_idx, inner_idx) = segment_index(segment, position)
    return unpack(segment.words[outer_idx], inner_idx, segment.granularity)
end

#Constructs the key corresoinding to the symbol
function cipher(symbol::Char, alphabet)::UInt
    i = 0
    for alphabet_symbol in alphabet
        if symbol == alphabet_symbol
            return UInt(i)
        end
        i += 1
    end
    error("Unrecognized symbol")
end

cipher_basis(symbol) = cipher(symbol, bases_symbols)
cipher_dotparen(symbol) = cipher(symbol, dotparen_symbols)

#Construct the symbol corresponding to the key
function decipher(key::UInt, alphabet)::Char
    return alphabet[key+1]
end

decipher_basis(key) = decipher(key, bases_keys)
decipher_dotparen(key) = decipher(key, dotparen_keys)

function print_segment(segment::Segment)
    println(String(decipher_basis.([unpack(segment, UInt(i)) for i in 1:segment.length])))
end

function invert(segment::Segment)::Segment
    invertd_segment = ZeroSegment(segment.length, segment.granularity)
    for i in 1:segment.length
        pack!(invertd_segment, unpack(segment, i), segment.length-i+1)
    end
    return invertd_segment
end

function complementary_segment(segment::Segment)::Segment
    Segment([~word for word in segment.words], segment.length, segment.granularity)
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

#SEQUENCES (Perfect Double Strands) ------------------------------------------------

function Sequence(string::String)::Sequence
    len = length(string)
    segment::Segment = ZeroSegment(len, 2)
    for i in 1:len
        pack!(segment, cipher_basis(string[i]), i)
    end
    return Sequence(segment)
end

function print_sequence(sequence::Sequence)
    print_segment(sequence.binary)
end

function invert(sequence::Sequence)::Sequence
    return Sequence(invert(sequence.binary))
end

function complementary_sequence(sequence::Sequence)::Sequence
    return Sequence(complementary_segment(sequence.binary))
end

function bitwise_binary_check(f::Function, sequence1::Sequence, sequence2::Sequence)::Bool
    return bitwise_binary_check(f, sequence1.binary, sequence2.binary)
end

are_equal(s1, s2) = bitwise_binary_check((a,b)->(a ⊻ b), s1, s2)

are_complementary(s1, s2) = bitwise_binary_check((a,b)->~(a ⊻ b), s1, invert(s2))

function inner_potentials(sequence::Sequence, data::DataFrame)::Vector{Float64}
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

function terminal_AT_penalty(sequence::Sequence)::Vector{Float64}
    segment = sequence.binary
    len = segment.length
    term_key = unpack(segment, len)
    if term_key == 0
        return [2.2e3, 6.9]
    end
    return [0.0, 0.0]
end

function symmetry_correction(sequence::Sequence)::Vector{Float64}
    if are_complementary(sequence, sequence)
        return [0.0, -1.4]
    end
    return [0.0, 0.0]
end

function sequence_potentials(sequence::Sequence, data_inner::DataFrame)::Vector{Float64}
    Δ = [0.2e3, -5.7]
    Δ += inner_potentials(sequence, data_inner)
    Δ += terminal_AT_penalty(sequence)
    Δ += symmetry_correction(sequence)
    return Δ
end

function melting_temperature(H, S, C)
    return H/(S+R*log(C/4))
end

function melting_curve(β, H, S)::Float64
    x = exp((β*H - S)/R)
    return 1 + x - sqrt(x^2 + 2x)
end

function create_NN_template(name::String)
    df = DataFrame(
        term = [String(decipher_basis.([unpack2(UInt(i), UInt(1)), unpack2(UInt(i), UInt(2))])) for i in 0:15],
        deltaH = [Float64(0) for _ in 0:15],
        deltaS = [Float64(0) for _ in 0:15],
    )
    CSV.write(name, df)
end

#DOUBLE STRANDS -----------------------------------------------------------------------------------

function pair(segment1::Segment, segment2::Segment)::Segment
    if segment1.length != segment2.length
        error("Segments must have same length to be paired")
    end
    length = segment1.length
    granularity = segment1.granularity + segment2.granularity
    inverted2 = invert(segment2)
    pair = ZeroSegment(length, granularity)
    for i in 1:length
        key = unpack(segment1, i) | (unpack(inverted2, i) << segment1.granularity)
        pack(pair, key, i)
    end
    return pair
end

function pair(sequence1::Sequence, sequence2::Sequence)::Sequence
    return Sequence(pair(sequence1.binary, sequence2.binary))
end

function Segment(key, granularity)
    return pack(ZeroSegment(1, granularity), key, 1)
end

function Sequence(key, granularity)
    return Sequence(Segment(key, granularity))
end

function subsegment(segment::Segment, idxs::Vector{T}) where T <: Integer
    len = length(idxs)
    subsegment = ZeroSegment(len, segment.granularity)
    for i in 1:len
        subsegment = pack(subsegment, unpack(segment, idxs[i]), i)
    end
    return subsegment
end

function subsequence(sequence::Sequence, idxs::Vector{T}) where T <: Integer
    return Sequence(subsegment(sequence.binary, idxs))
end

function extend(segment::Segment, key::UInt)::Segment
    len = segment.length + 1
    extended_segment = ZeroSegment(len, segment.granularity)
    for i in 1:segment.length
        pack(extended_segment, unpack(segment, i), i)
    end
    pack(extended_segment, key, len)
    return extended_segment
end

function extend(sequence::Sequence, key::UInt)::Sequence
    return Sequence(extend(sequence.binary, key))
end


function extend!(hairPin::HairPin, pair::UInt)::HairPin
    hairPin.paired_sequence.binary = extend(hairPin.paired_sequence.binary, pair)
    return hairPin
end

#INTERACTING SINGLE STRANDS -----------------------------------------------------------------------

function substrand(strand::SingleStrand, idxs::Vector{T}) where T <: Integer
    return SingleStrand(Sequence(subsegment(strand.sequence.binary, idxs)), Structure(subsegment(strand.structure.binary, idxs)), length(idxs))
end

function SingleStrand(sequence_lit::String, dotparen_lit::String)
    len = length(sequence_lit)
    if len != length(dotparen_lit)
        error("Sequences and secondary structures must have the same length")
    end
    sequence_binary::Segment = ZeroSegment(len, 2)
    structure_binary::Segment = ZeroSegment(len, 2)
    for i in 1:len
        pack!(sequence_binary, cipher_basis(sequence_lit[i]), i)
        pack!(structure_binary, cipher_dotparen(dotparen_lit[i]), i)
    end
    return SingleStrand(Sequence(sequence_binary), Structure(structure_binary), len)
end

function Ring(length, pair::UInt)
    HairPin(length, Sequence(pair, 4))
end

function integrate(structure::Structure)
    integral = Vector{UInt}(undef, structure.binary.length+1)
    integral[1] = UInt(0)
    for i in 1:structure.binary.length
        integral[i+1] = integral[i] + unpack(structure.binary, i) - 1
    end
    return integral
end 


function extract_hairpins!(strand::SingleStrand, hairPins::Vector{HairPin})
    integral = integrate(strand.structure)
    if integral[end] != 0
        error("Invalid structure")
    end
    ones_idxs = findall(x -> x == 1, integral)
    loop_length = count(x -> x == UInt(1), [unpack(strand.structure.binary, UInt(i)) for i in ones_idxs])
    first_element = unpack(strand.sequence.binary, 1)
    second_element = unpack(strand.sequence.binary, strand.length)
    println(bitstring(first_element))
    println(bitstring(second_element))
    println()
    if ~(first_element ⊻ second_element) & masks[1] == 0
        outer_pair = first_element | second_element << 2
    else
        error("Structure and sequence are not compatible")
    end
    attachment_number = 0
    i = 1
    while i <= length(ones_idxs)
        structure = unpack(strand.structure.binary, UInt(ones_idxs[i]))
        if structure == 2
            extract_hairpins!(substrand(strand, collect(ones_idxs[i]:ones_idxs[i+1]-1)), hairPins)
            attachment_number += 1
        end
        i += 1
    end
    println(attachment_number)
    if loop_length == 0 && attachment_number == 1
        extend!(hairPins[end], outer_pair)
    else 
        push!(hairPins, Ring(loop_length, outer_pair))
    end
end

function loop_contribution_G(hairPin::HairPin)
    if 2 < hairPin.loop_length < 10
        loop_contribution_array = [3.2,3.6,4.0,4.4,4.6,4.7,4.8]
        return loop_contribution_array[hairPin.loop_length]
    elseif hairPin.loop_length > 9
        return 4.8 + 1.75*R*T
    else
    end 
end

function loop_contribution_H()
end

function loop_contribution_S(T::Float64)
end

function hairpin_potentials(hairPin::HairPin, data_inner::DataFrame)::Vector{Float64}
    loop_contributions_G37 = [3.2,3.6,4.0,4.4,4.6,4.7,4.8]
    loop_contributions_H = []
    Δ = sequence_potentials(hairPin.paired_sequence, data_inner)
    
end

strand = SingleStrand("TA","()")

hairPins = Vector{HairPin}([])

extract_hairpins!(strand, hairPins)

hairPins



end