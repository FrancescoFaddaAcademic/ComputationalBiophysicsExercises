using HTTP
using StaticArrays
using CairoMakie

const FS = Sys.WORD_SIZE #Frame size (size of the partition blocks, must consider system architecture)
const EPF = Int(FS/2) #Number of bases in a frame 

struct Strand #A strand always read from 5' to 3'
    length::UInt
    bytes::UInt
    literal::String
    binary::Vector{UInt64}
end

#Matched pairs potentials

#a g c t/u

const R = 1.987 #cal/K/mol
const C = 1

const ΔH = Vector{Float64}([-7.6, -8.2, -8.5, -7.2, -7.8, -8.0, -10.6, -8.5, -8.4, -9.8, -8.0, -8.2, -7.2, -8.4, -7.8, -7.6] * 1000)
const ΔS = Vector{Float64}([-21.3, -22.2, -22.7, -21.3, -21.0, -19.9, -27.2, -22.7, -22.4, -24.4, -19.9, -22.2, -20.4, -22.4, -21.0, -21.3])

const H_init::Float64 = 0.2 
const S_init::Float64 = -5.7

const H_TermATpen::Float64 = 2.2
const S_TermATpen::Float64 = 6.9

const H_symCor::Float64 = 0.0
const S_symCor::Float64 = -1.4

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

key00 = 0x0000000000000000
key01 = 0x5555555555555555
key10 = 0xaaaaaaaaaaaaaaaa
key11 = 0xffffffffffffffff

Strand(string::String) = Strand(length(string), UInt64(ceil(length(string)/EPF)), string, litToBin(string)[2]) 

reverseStrand(strand::Strand) = Strand(reverse(strand.literal))

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

function areComplementary(strand1::Strand, strand2::Strand)::Bool
    if strand1.length ≠ strand2.length
        return false
    end
    nBytes = strand1.bytes 
    for i in 1:(nBytes-1)
        if strand1.binary[i] ⊻ reverseStrand(strand2).binary[i] ≠ ~0x0000000000000000
            return false
        end
    end
    if (strand1.binary[nBytes] ⊻ reverseStrand(strand2).binary[nBytes] ⊻ sumMasks[strand1.length - (nBytes-1)*EPF]) ≠ 0x0000000000000000
        return false
    end
    return true
end

function bindingThermoComp(strand::Strand)::SVector{2,Float64}

    #Initiation
    H = 0.2
    S = -5.7

    #Evaluation of the inner terms
    for oId in 1:strand.bytes-1
        for iId in 1:EPF-1
            bind_type_Id = (strand.binary[oId] & (masks[iId]|masks[iId+1])) >> (2*(iId-1))
            H += ΔH[bind_type_Id+1]
            S += ΔS[bind_type_Id+1]
        end
        bind_type_Id = (strand.binary[oId] & masks[EPF] >> (2*(EPF-1))) | (strand.binary[oId + 1] & masks[1] << 2)
        H += ΔH[bind_type_Id+1]
        S += ΔS[bind_type_Id+1]
    end
    for iId in 1:(strand.length-(strand.bytes-1)*EPF-1)
        bind_type_Id = (strand.binary[strand.bytes] & (masks[iId]|masks[iId+1])) >> (2*(iId-1))
        H += ΔH[bind_type_Id+1]
        S += ΔS[bind_type_Id+1]
    end

    #Evaluation of the border terms
    if strand.literal[1] ∈ ['a', 'A', 't', 'T']
        H += 2.2
        S += 6.9
    end
    if strand.literal[strand.length] ∈ ['a', 'A', 't', 'T']
        H += 2.2
        S += 6.9
    end

    #Symmetry Correction
    if areComplementary(strand, strand)
        S += -1.4
    end

    @show return SVector{2,Float64}([H,S])
end

function meltingTemperature(H, S, C)
    return H/(S+R*log(C/4))
end

function meltingCurve(β, H, S)::Float64
    x = exp(-(β*H - S)/R)
    return 1-2(-1+sqrt(1+2*x))/(2*x)
end

println("-----------------------------------------------------------")

#CGTTGA


strand = Strand("CGTTGA")

data = bindingThermoComp(strand)

H = data[1]
S = data[2]
Tm = meltingTemperature(data...,C)

print("Strand = ", strand.literal, "\n")

print("ΔH = ")
print(round(H, digits = 1))
print(" cal/mol \n")

print("ΔS = ")
print(round(S, digits = 1))
print(" cal/mol \n")

print("Tm = ")
print(Tm-273.15)
print(" C° (")
print(Tm)
print(" K)\n")

Ts= Vector(0:1:600)
Xs = Vector{Float64}(undef, length(Ts))

for i in 1:length(Xs)
    Xs[i] = meltingCurve(1/Ts[i], H, S)
end

meltingCurve(1/Tm,H,S)

lines(Ts, Xs)

println("-----------------------------------------------------------")