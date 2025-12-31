# SantaLucia.jl

### Implemented extensions 
- Variable Strand Concentration
- Melting curve plot
- Internal loops, bulges and Hairpins (Internal loops and Hairpins must however be longer than 4 and 3 respectively due to the complexity of the data for the shortest case)

## Documentation

### Some preliminary notes on the Bitwise Manipulation
To keep memory and computations efficient I've decided to use a bitwise data control, in hindsight this contributes to a worse legibility of the code, however I think that the efficiency advantages together with the improved elegance are worth the slightly more complex syntax. 

In the following I will breefly discuss some of the methods that I've implemented for the bit manipulation
```julia
mutable struct Segment
    words::Vector{UInt}
    length::UInt
    granularity::UInt
end
```
This is the main structure that I've used all over the code to store information both for the primary and secondary structures. A `Segment` is a sequence of `Words`, a word is nothing but a raw piece of memory with a specified size, in most modern computer such size is 64 bits. A segment serves as a very customizable memory structure, it will contain `elements` of some length specified by the `granularity`, this granularity must be smaller than the size of a single word (at least by a factor of two to make the structure even remotely useful). For example in the case of sequences of nucleic acids each element is one of four letter and hence the granularity of this kind of data will be 2 (2^2 = 4). For secondary structure there are only 3 symbols ('(', '.', ')'), however one must chose 2 as the granularity also in this case as it is the smallest power of two that enables the storage of each symbol. This means that for sequences of length <= 32 a single word for the primary structure and a single word for the secondary structure. Lastly the lenght is just the amount of elements contained in the segmets.

The reason this approach is so appropriate for this usecase is that most operations are done equally on all element of a structure, this means that through bitwise manipulation one is able to "parallelise" the computation and reduce the memory allocated up to a factor of 32! While this is irrelevant for the size of the sequences processed in this assignment, however for longer sequences it car represent a good computational advantage. 

### Segment index
```julia
function segment_index(segment::Segment, position::Integer)::Tuple{UInt, UInt}
    word_capacity = div(word_size, segment.granularity)
    outer_idx = cld(position, word_capacity)
    inner_idx = mod1(position, word_capacity)
    return (outer_idx, inner_idx)
end

```
To make all of this work one must implement some methods for accessing the elements stored in the various segments. Let's start by the `segment_index` function, this is a simple funciton that takes as an input a segment and a position and returns a tuple corresponding to the word index and the inner position of that index.

### Pack and Unpack
```julia
function pack(word::UInt, key::UInt, position::Integer, granularity::Integer)::UInt
    return key << (granularity*(position-1)) | (word & ~masks[position])
end

function unpack(word::UInt, position::Integer, keysize::Integer)::UInt
    return (word & masks[position]) >> (keysize*(position-1))
end
```
These functions are used respectively to put and read elements (`key`) into words. There are analogous functions for segments that use the previous `segment_index` function to select the correct word and position before using these to do the actual data manipulation. In `pack` the key is bitshifted to the desired position, then substituted to the original bits (a mask is nothing but a word with zeros everywhere but in the desired position). Similarly in `unpack` one first selects the bit with the appropriate mask and then bitshifts them to the beginning to extract the key.

### Cipher and Decipher
```julia
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

function decipher(key::UInt, alphabet)::Char
    return alphabet[key+1]
end
```
To interface beteen the human readable data and the machine efficient data some interface must be set up, the functions `cipher` and `decipher` simply put in biunivoque correspondence elements of our alphabets (e.g. 'A','G','C','T' for the primary structure) to elements (or keys). This correspondence is often done in a specific way to implement the following functionality in the most natural manner. For example, regarding the primary structure the correspondence will be the following: 
- A -> 00, 
- G -> 01,
- C -> 10,
- T -> 11

By doing this one can check if two bases are complementary by simply performing a XOR and then checking if they are 0, as I was saying this kinds of operations can be performed on each word simultaneously.
In the code there are other 

### Core structures
#### Comptime constants
```julia 
const word_size = Sys.WORD_SIZE
const bases_symbols = SVector{4,Char}(['A', 'G', 'C', 'T'])
const dotparen_symbols = SVector{3, Char}([')','.','('])
const keys = SVector{4,UInt}([0, 1, 2, 3])
const masks = SVector{32,UInt}([2^(2*i)+2^(2*i + 1) for i in 0:(word_size/2)-1])
const R = 1.987
```
These are some useful constants that are used all over the code

#### Sequence
```julia
mutable struct Sequence
    binary::Segment
end
```
#### Single Strand
```julia
mutable struct SingleStrand
    sequence::Sequence
    structure::Structure
    length::UInt
end
```

#### Pin
```julia
mutable struct Pin
    loop_length::UInt
    paired_sequence::Sequence
    type::Char
end
```

### Core functions

#### Complementarity of Sequences
```julia
are_complementary(s1, s2) = bitwise_binary_check((a,b)->~(a ⊻ b), s1, invert(s2))

```
#### Inner Thermodynamic Contribution
```julia 
function inner_thermodynamic_contribution(sequence::Sequence, data::DataFrame)::Vector{Float64}
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
```
#### Terminal AT Penalty
```julia
function terminal_AT_penalty(sequence::Sequence)::Vector{Float64}
    segment = sequence.binary
    len = segment.length
    term_key = unpack(segment, len)
    if term_key == 0
        return [2.2e3, 6.9]
    end
    return [0.0, 0.0]
end
```

#### Symmetry Correction
```julia
function symmetry_correction(sequence::Sequence)::Vector{Float64}
    if are_complementary(sequence, sequence)
        return [0.0, -1.4]
    end
    return [0.0, 0.0]
end

```
#### Sequence Thermodynamics
```julia
function sequence_potentials(sequence::Sequence, data_inner::DataFrame)::Vector{Float64}
    Δ = [0.2e3, -5.7]
    Δ += inner_thermodynamic_contribution(sequence, data_inner)
    Δ += terminal_AT_penalty(sequence)
    Δ += symmetry_correction(sequence)
    return Δ
end
```

#### Melting temperature 
```julia
function melting_temperature(H, S, C)
    return H/(S+R*log(C/4))
end
```

#### Melting curve
```julia
function melting_curve(β, H, S)::Float64
    x = exp((β*H - S)/R)
    return 1 + x - sqrt(x^2 + 2x)
end
```

#### Pin Extraction
```julia
function extract_pins!(strand::SingleStrand, pins::Vector{Pin})
    integral = integrate(strand.structure)
    if integral[end] != 0
        error("Invalid structure")
    end
    ones_idxs = findall(x -> x == 1, integral)
    loop_length = count(x -> x == UInt(1), [unpack(strand.structure.binary, UInt(i)) for i in ones_idxs])
    first_element = unpack(strand.sequence.binary, 1)
    second_element = unpack(strand.sequence.binary, strand.length)
    if ~(first_element ⊻ second_element) & masks[1] == 0
        outer_pair = first_element | second_element << 2
    else
        error("Structure and sequence are not compatible")
    end
    has_marginal = false
    children_number = 0
    i = 1
    while i <= length(ones_idxs)
        structure = unpack(strand.structure.binary, ones_idxs[i])
        if structure == 2
            if unpack(strand.structure.binary, ones_idxs[i]-1) == 2 || unpack(strand.structure.binary, ones_idxs[i+1]) == 0
                has_marginal = true
            end
            extract_pins!(substrand(strand, collect(ones_idxs[i]:ones_idxs[i+1]-1)), pins)
            children_number += 1
        end
        i += 1
    end
    only_child = children_number == 1

    if loop_length == 0 && only_child 
        extend!(pins[end], outer_pair)
    elseif has_marginal && only_child
        push!(pins, Ring(loop_length, outer_pair, 'B'))
    elseif children_number == 0 
        push!(pins, Ring(loop_length, outer_pair, 'H'))
    else
        push!(pins, Ring(loop_length, outer_pair, 'I'))
    end
end
```

#### Loop Thermodynamic Contribution
```julia
function loop_thermodynamic_contribution(pin::Pin, hairpin_data::DataFrame, internal_loop_data::DataFrame, bulge_data::DataFrame)::Vector{Float64}
    if pin.type == 'H'
        is_long = pin.loop_length > 9
        data = hairpin_data
        if pin.loop_length < 3
            error("Hairpin of length < 3 are not supported")
        elseif is_long
            return [data[9,"Enthalpy"], data[9,"Entropy"]+1.75*R*log(pin.loop_length/9)]
        else
            return [data[pin.loop_length,"Enthalpy"], data[pin.loop_length,"Entropy"]]
        end
    elseif pin.type == 'I'
        is_long = pin.loop_length > 6
        data = internal_loop_data
        if pin.loop_length < 4
            error("Internal loop of length < 4 are not supported")
        elseif is_long
            return [data[6,"Enthalpy"], data[9,"Entropy"]+1.08*R*log(pin.loop_length/6)]
        else
            return [data[pin.loop_length,"Enthalpy"], data[pin.loop_length,"Entropy"]]
        end
    elseif pin.type == 'B'
        is_long = pin.loop_length > 9
        data = bulge_data
        if is_long
            return [data[6,"Enthalpy"], data[9,"Entropy"]+1.75*R*log(pin.loop_length/6)]
        else
            return [data[pin.loop_length,"Enthalpy"], data[pin.loop_length,"Entropy"]]
        end
    else
        error("Unrecognized loop structure")
    end
    if is_long
        return [data[9,"Enthalpy"], data[9,"Entropy"]+1.75*R*log(pin.loop_length/9)]
    else
        return [data[pin.loop_length,"Enthalpy"], data[pin.loop_length,"Entropy"]]
    end
end
```

## Additional notes

