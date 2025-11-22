using HTTP
using StaticArrays

struct Strand #A strand always read from 5' to 3'
    length::UInt32
    literal::String
    binary::Vector{UInt8}
end

global masks = SVector{4,UInt8}(3, 12, 48, 192)
global sumMasks = SVector{4,UInt8}(3, 15, 63, 255)

Strand(string::String) = Strand(length(string), string, litToBin(string)[2]) 

reverse(strand::Strand) = Strand(reverse(strand.literal))

function litToBin(lit::String)
    len = length(lit)
    nBytes = UInt8(ceil(len/4))
    bins = zeros(UInt8,nBytes)
    masks = SVector{4,UInt8}(3, 12, 48, 192)
    isR::Bool = false
    isD::Bool = false
    for oId in 1:nBytes #Outer loop on the bytes
        for iId in 1:4       #Inner loop on the bases 
            litId = (oId-1)*4+iId
            if litId > len; break end #Bound check
            char = lit[litId]
            mask = masks[iId]
            #Purines
            if char == 'a' || char == 'A'
                key = 0
            elseif char == 'g' || char == 'G'
                key = 85
            #Pyrimidines
            elseif char == 'c' || char == 'C'
                key = 170    
            elseif char == 'u' || char == 'U'
                key = 255
                isR = true
            elseif char == 't' || char == 'T'
                key = 255
                isD = true
            else
                error("Invalid Literal Sequence")
            end
        bins[oId] |= key & mask
        end
    end
    if isD && isR
        error("Invalid Literal Sequence: Both Uracil and Thyamine present")
    else
        return (len,bins)
    end
end

function areComplementary(strand1::Strand, strand2::Strand)::Bool
    if strand1.length != strand2.length
        return false
    end
    nBytes = UInt32(ceil(strand1.length/4))
    for i in 1:(nBytes-1)
        if strand1.binary[i] ⊻ reverse(strand2).binary[i] != !0
            return false
        end
    end
    if (strand1.binary[nBytes] ⊻ reverse(strand2).binary[nBytes] ⊻ sumMasks[strand1.length - (nBytes-1)*4]) != 0
        return false
    end
    return true
end

print(bitstring.(litToBin("agct")[2]))

areComplementary(Strand("agctactg"), Strand("agctactg"))