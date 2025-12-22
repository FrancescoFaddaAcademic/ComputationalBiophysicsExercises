using StaticArrays

const FS = 64 #Frame size (size of the partition blocks, must consider system architecture)
const EPF = Int(FS/2) #Number of bases in a frame 

struct Sec_struct
    dotparen::String
    dot::Vector
    open::Vector
    closed::Vector
end

struct Strand
    length::UInt64
    bytes::UInt64
    literal::String
    binary::Vector{UInt64}
    sec_struct::Sec_struct
end

prov_masks = Vector{UInt64}(undef, EPF) 
prov_sum_masks = Vector{UInt64}(undef, EPF)
prov_masks[1] = prov_sum_masks[1] = 0x0000000000000003

for i in 2:EPF
    local mask = UInt64(0x0000000000000003 << (2*(i-1)))
    prov_masks[i] = mask
    prov_sum_masks[i] = prov_sum_masks[i-1] + mask
end

const masks = SVector{EPF,UInt64}(prov_masks)
const sumMasks = SVector{EPF,UInt64}(prov_sum_masks)

const key00 = 0x0000000000000000
const key01 = 0x5555555555555555
const key10 = 0xaaaaaaaaaaaaaaaa
const key11 = 0xffffffffffffffff

function Strand(string::String, dotparen::String)
    if length(string) == length(dotparen) 
        Strand(length(string), UInt64(ceil(length(string)/EPF)), string, litToBin(string)[2], dotparen_parser(dotparen)) 
    else
        error("The sequence literal and structure proposed have different length, hence, they are incompatible")
    end
end
reverseStrand(strand::Strand, dotparen::String) = Strand(reverse(strand.literal), dotparen)

function litToBin(lit::String)
    len = length(lit)
    nBytes = UInt8(ceil(len/EPF))
    bins = zeros(UInt64,nBytes)
    isR::Bool = false
    isD::Bool = false
    for oId in 1:nBytes #Outer loop on the bytes
        for iId in 1:EPF       #Inner loop on the bases 
            litId = (oId-1)*EPF+iId
            if litId > len; break end #Bound check
            char = lit[litId]
            mask = masks[iId]
            #Purines
            if char ∈ ['a','A']
                key = key00
            elseif char ∈ ['g','G']
                key = key01
            #Pyrimidines
            elseif char ∈ ['c','C']
                key = key10    
            elseif char ∈ ['u','U']
                key = key11
                isR = true
            elseif char ∈ ['t','T']
                key = key11
                isD = true
            else
                error("Invalid Literal Sequence")
            end
        bins[oId] |= key & mask
        end
    end
    if isD && isR
        error("Invalid sequence: Both Uracil and Thyamine present")
    end
    return (len,bins)
end

function dotparen_parser(dotparen::String)::Sec_struct
    len = length(dotparen)
    prov_dot_idxs = Vector{UInt64}(undef, len)
    prov_open_idxs = Vector{UInt64}(undef, len)
    prov_closed_idxs = Vector{UInt64}(undef, len)
    dot_number::UInt64 = 0
    open_number::UInt64 = 0
    closed_number::UInt = 0
    for idx in 1:len
        char = dotparen[idx]
        if char == '.'
            dot_number += 1
            prov_dot_idxs[dot_number] = idx
        elseif char == '('
            open_number += 1
            prov_open_idxs[open_number] = idx
        elseif char == ')'
            closed_number += 1
            prov_closed_idxs[closed_number] = idx
        else
            error("Improper content: Only '.' '(' ')' characters are admitted")
        end
    end
    if open_number ≠ closed_number
        error("Improper content: different number of open and closed parentheses")
    end
    if dot_number ≠ 0  
        dot_idxs = prov_dot_idxs[1:dot_number]
    else
        dot_idxs = Vector{UInt64}(undef, 0)
    end
    if open_number ≠ 0
        open_idxs = prov_open_idxs[1:open_number]
    else
        open_idxs = Vector{UInt64}(undef, 0)
    end
    if closed_number ≠ 0
        closed_idxs = prov_closed_idxs[1:closed_number]
    else
        closed_idxs = Vector{UInt64}(undef, 0)
    end
    return Sec_struct(dotparen, dot_idxs, open_idxs, closed_idxs)
end

function check_compatibility(sec_struct::Sec_struct)
    bound_len = length(sec_struct.closed)
    lit_pairs = Vector{SVector{2, UInt64}}(undef, bound_len)
    pair_idx = 1
    if bound_len == 0
        return true
    end
    first_closed = sec_struct.closed[1]
    for open_idx in 1:bound_len
        if sec_struct_open[open_idx] > first_closed
            if open_idx == 1
                error("The provided structure is invalid")
            else
                lit_pairs[lit_idx] = SVector([first_closed, sec_struct_open[open_idx-1]])
                unpaired_len = sec_struct_open[open_idx-1] - first_closed - 1 

            end 
        end
    end
end

function remove_window(s::String, min::T, max::S)::String where {T, S<:Integer}
    s[1:min-1] * s[max+1:length(s)]
end

function recursion!(orig_strand::Strand, curr_strand::Strand, paired_sections::Vector{Vector{SVector{2,Char}}}, loops::Vector{UInt64}; iteration = 1, sec_idx = 1, pair_idx = 1, print = true)
    bound_len = length(curr_strand.sec_struct.closed)
    if bound_len == 0
        return
    end
    first_closed_idx = curr_strand.sec_struct.closed[1]
    if curr_strand.sec_struct.open[1] > first_closed_idx
        error("Structure error: closed parenthesis before any open")
    end
    matched_open_idx = curr_strand.sec_struct.open[bound_len]
    if bound_len == 1
        matched_open_idx = curr_strand.sec_struct.open[1]
    else
        for i in 1:bound_len-1
            next_open_idx = curr_strand.sec_struct.open[i+1]
            if next_open_idx > first_closed_idx
                matched_open_idx = curr_strand.sec_struct.open[i]
                break
            end
        end
    end
    loop_len = first_closed_idx - matched_open_idx - 1
    loops[orig_strand.sec_struct.closed[iteration]] = loop_len
    if iteration == 1
    elseif loop_len ≠ 0
        sec_idx += 1
        pair_idx = 1
    else
        pair_idx += 1
    end
    paired_sections[sec_idx][pair_idx] = SVector{2,Char}([curr_strand.literal[matched_open_idx],curr_strand.literal[first_closed_idx]])
    new_literal = remove_window(curr_strand.literal, matched_open_idx, first_closed_idx)
    new_sec_struct = remove_window(curr_strand.sec_struct.dotparen, matched_open_idx, first_closed_idx)
    new_curr_strand = Strand(new_literal, new_sec_struct)
    if print
        println(curr_strand.sec_struct.dotparen* " --> loop length = $loop_len, pair = $(curr_strand.literal[matched_open_idx]) - $(curr_strand.literal[first_closed_idx])")
    end
    recursion!(orig_strand, new_curr_strand, paired_sections, loops; iteration = iteration + 1, sec_idx = sec_idx, pair_idx = pair_idx)
end

dotparen = "((((.).).(...)))"
literal = "aaaaaaaaaaaaaaaa" 

len = length(dotparen)

strand = Strand(literal, dotparen)

paired_sections = Vector{Vector{SVector{2,Char}}}(undef, len)
for i in 1:len
    paired_sections[i] = fill(SVector{2,Char}(['.','.']), len)
end

loops = fill(UInt64(0), len)

recursion!(strand, strand, paired_sections, loops)

#Note: La parte di parsing dotparen può essere parallelizzata come albero 
#se si tiene conto della profondità del dominio onde evitare i problemi di concurrency.