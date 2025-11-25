using StaticArrays

struct Sec_struct
    dotparen::String
    dot::Vector
    open::Vector
    closed::Vector
end

dotparen = "(((...((...)))))"

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
        closed_idxs = prov_closed_idxs[1:dot_number]
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

function recursive_reduction(sec_structure::Sec_struct, lit_pairs::Vector{SVector{2, UInt64}}, loops::Vector{UInt64}, pair_idx::Int, loop_idx::Int)::Sec_struct
    bound_len = length(sec_struct.closed)
    pair_idx = 1
    if bound_len == 0
        return
    end
    first_closed = sec_struct.closed[1]
    for i in 1:bound_len
        if sec_struct_open[i] > first_closed
            if open_idx == 1
                error("The provided structure is invalid")
            else
                lit_pairs[lit_idx] = SVector([first_closed, sec_struct_open[i-1]])
                loops[loop_idx] = sec_struct_open[i-1] - first_closed - 1 
                lit_idx += 1
                loop_idx += 1
                #Here I have to implement the recursion
            end 
        end
    end
end