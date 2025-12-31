# SantaLucia
## Brief theoretical introduction and architecture of the code
The SantaLucia model is the most widely used Nearest Neighbor model for dna, it is used to calculate the thermodynamic potentials associated with the formation of the secondary structure of some nucleic acid. The main idea consists in the assumption that the most relevant interaction is present at the level of adjacent pairs, each giving an intependent contribution to the thermodynamic potentials of the whole system by this assumption of separability one manages to achieve low computational cost while still preserving most of the relevant thermodynamic information. In general a complete NN model can be quite complex, taking into account many of the possible structure motifs that may arise from interaction. The main two kinds of contributions that such an NN model calculates are the following
- The adjacent nucleotide-pair interaction (in a section that is paired), together with some corrections
- The contribution of unpaired sections giving rise to more complex geometries (such as Bulges, Hairpins, Internal-Loops etc.)  

In the present assignment I've decided to model only a subclass of these motifs, in particular every paired section is assumed to be complete (the code will throw an error if you try to match two non complementary strands), and while Bulges Hairpins and Interal-Loops are supported for the latter two there is a minimum size needed for the code to work, in particular a hairpin should have a size of at least three bases while the internal loop of at least 4. I've avoided these particular cases to keep the code and (mostly) data realitively simple, as, due to the close configuration there is no easy rule for these cases, and a lot of unique data is needed for a correct modeling. In short the implemented extensions are then the following  

### Implemented extensions 
- Variable Strand Concentration
- Melting curve plot
- Internal loops, bulges and Hairpins (Internal loops and Hairpins must however be longer than 4 and 3 respectively due to the complexity of the data for the shortest case)

### The Mathematics of RNA Secondary structure
RNA secondary structure can be quite complex, to understand it, and to correctly compute its thermodynamics, i've found that it is useful to look at any RNA secondary structure as a sort of composable graph. To discuss this, let's start by looking at the dotparen notation. The dotparen notation used to concisely describe the secondary structure of a piece of RNA, each basis is associated to one of the following symbols: '(', '.' or ')'. Opposing parenthesis represent paired nucleotides in the structure, while dots unpaired ones. While this is sufficient to completely describe the secondary structure (at least at the degrees of freedom with higher energy), however.  

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

By doing this one can check if two bases are complementary by simply performing a XNOR and then checking if they are 0, as I was saying this kinds of operations can be performed on each word simultaneously.
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
These are some useful constants that are used all over the code, mostly to do the bit manipulation.

#### Sequence
```julia
mutable struct Sequence
    binary::Segment
end
```
At the moment a sequence coincides with its binary sequence, I've decided to create this encapsulation to be able to add some additional information if needed. 
NOTE: Due to the lack of a support for mismatches every "perfect" double sequence will be treated just as a single sequence as including the second sequence would have been redundant. 

#### Single Strand
```julia
mutable struct SingleStrand
    sequence::Sequence
    structure::Structure
    length::UInt
end
```
A single strand contains both its sequence, its structure (that is analogous to the sequence structure) and its length

#### Pin
```julia
mutable struct Pin
    loop_length::UInt
    paired_sequence::Sequence
    type::Char
end
```
In this instance, where the interaction is dependent only on the size and on the type of loop, a Pin is well described by the length and type of its loop, and by the paired sequence (that is a sequence in which each element is a couple of connected bases).

### Core functions

#### Complementarity of Sequences
```julia
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

are_complementary(s1, s2) = bitwise_binary_check((a,b)->~(a ⊻ b), s1, invert(s2))

```
To check if two sequences are complementary we start by swapping the orientation of one of them using invert, then we use a specific instance of the `bitwise_binary_check` to parallelize the bitwise operation described by the lambda `(a,b)->~(a ⊻ b)` (that is simply an XNOR) on the entire content of the binary sequence. Due to the choice of "cifration" it can be checked that this operation returns 0 if and only if the two strings are complementary.
NOTE: To check complementarity one has to insert both sequences in the same orientation, for example 5' -> 3'
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
This function uses the cifration to access the thermodynamic data provided in the `data/` directory (for htis particular applliaction the `Watson-Crick-NN-Parameters.csv` file). The idea here is that each couple defines uniquely a number beteen 0 and 15. For example, if one takes AG to be a couple of nearest neighbors, this couple would appear as 0100 in the corresponding segment, this is 4 in binary. By listing enthalpy and ertropy values in a coherent fashion the elements can be used directly to access the dataframe at the correct location. It is important to note that there is a redundance as there is no difference looking at one strand or the other. In particular, still running 5' -> 3', opposite to the previous couple there will be CT, that corresponds to 1110 (that is 14) so these tho rows will be identical. To be able to support mismatches the data structures should be modified, in particular I believe the better option would be to create a struct containing a sequence whose elements are copuples of opposing bases. by doing so a similar approach could be used, this time using a table with 256 rows.
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
This function simply checks if the last element of the sequence is an 'A', and applies the thermodinamic modifications accordingly
#### Symmetry Correction
```julia
function symmetry_correction(sequence::Sequence)::Vector{Float64}
    if are_complementary(sequence, sequence)
        return [0.0, -1.4]
    end
    return [0.0, 0.0]
end

```
This function checks if the sequence is self complementary and applies the thermodynamic modifications accordingly
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
Lastly this function composes all the contribution of the various terms giving the final thermodynamic potentials for the double strand (and as we will see for some of the single strand).
#### Melting temperature and Melting Curve 
```julia
function melting_temperature(H, S, C)
    return H/(S+R*log(C/4))
end

function melting_curve(β, H, S)::Float64
    x = exp((β*H - S)/R)
    return 1 + x - sqrt(x^2 + 2x)
end
```
These two functions compute the melting temperature and melting curve, in the latter the form is different to that proposed to avoid precision errors due to the exponential function both at the numerator and denominator in the original version.

#### Pin Extraction
```julia
function extract_pins!(strand::SingleStrand, pins::Vector{Pin})
    integral = integrate(strand.structure)
    if integral[end] != 0
        error("Invalid structure")
    end
#1
    ones_idxs = findall(x -> x == 1, integral)
    loop_length = count(x -> x == UInt(1), [unpack(strand.structure.binary, UInt(i)) for i in ones_idxs])
    first_element = unpack(strand.sequence.binary, 1)
    second_element = unpack(strand.sequence.binary, strand.length)
#2
    if ~(first_element ⊻ second_element) & masks[1] == 0
        outer_pair = first_element | second_element << 2
    else
        error("Structure and sequence are not compatible")
    end
    has_marginal = false
    children_number = 0
    i = 1
#3
    while i <= length(ones_idxs)
        structure = unpack(strand.structure.binary, ones_idxs[i])
        if structure == 2
#4
            if unpack(strand.structure.binary, ones_idxs[i]-1) == 2 || unpack(strand.structure.binary, ones_idxs[i+1]) == 0
                has_marginal = true
            end
#5
            extract_pins!(substrand(strand, collect(ones_idxs[i]:ones_idxs[i+1]-1)), pins)
            children_number += 1
        end
        i += 1
    end
    only_child = children_number == 1
#6
    if loop_length == 0 && only_child 
        extend!(pins[end], outer_pair)
    elseif has_marginal && only_child
        push!(pins, Ring(loop_length, outer_pair, 'B'))
#7
    elseif children_number == 0 
        push!(pins, Ring(loop_length, outer_pair, 'H'))
    else
        push!(pins, Ring(loop_length, outer_pair, 'I'))
    end
end
```
This is the main decomposition function, it starts by taking the secondary structure as a segment, the cifration of the dotparen string is done as follows:
- ')' -> 00 = 1 + (-1)
- '.' -> 01 = 1 + (0)
- '(' -> 10 = 1 + (+1)

By subtracting 1 to each value one can think of the dotparen notation as a sort of derivative of a function, in particular one can antiderive to obtain a landscape that describes at what level each element is in the structure. In particular, paired bases, and elements of the same loop should be at the same level. 

At this point the function proceeds as follows (numbers are presented as comment in the code for readability): 
1. Looks at where the ones are (this will correspond to the level immediately deeper).
2. Stitches the pair corresponding to the outermost parentheses and checks that the pair is complementary. 
3. Runs through the position of the ones and assess what is on their right (a '.' or a '(') and how many of each there are.
4. If it is a '(' it checks if it (or the next) is the first (or the last) one, if so it enables a boolean flag `has_marginal`
5. If it is a '(' it calles itself on the higher level recursively.
6. Decides what kind of pin it is and assigns the correspondin identifier as a `type`.
7. If there are no points, insead of creating another Pin, the last one is extended.

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
Finally this function extracts the correct values for the enthalpy and entropy associated to the loop type and lenght from the corresponding databases or calculates them for longer loop lenght.
## Additional notes



