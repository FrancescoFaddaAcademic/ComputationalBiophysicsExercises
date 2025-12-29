using HPModel
using LaTeXStrings
using GLMakie
using BenchmarkTools
using Profile

#This test reproduces Table1 and Figure2 from the article [Chan-Dill 1989] 

function plot_compactness_categorized(categorized_paths::Vector{Vector{Path}})
    bins = Base.length(categorized_paths)
    binwidth = 1/(bins-1)
    fig = Figure(size=(800, 600), fontsize = 25)
    ax = Axis(
        fig[1, 1], 
        #title = "Number of sequences by compactness", 
        xlabel = L"\rho", 
        ylabel = "Number of sequences",
        xgridvisible=false,
        ygridvisible=true,
        xticklabelsize = 20,
        yticklabelsize = 18,
        xticks = round.(0:1/(bins-1):1, digits = 2),
        yticklabelrotation = π/4,
        )
    xlims!(ax, -binwidth/2, 1+binwidth/2)
    ylims!(ax, 0, nothing)
    barplot!(ax, (0:bins-1)./(bins-1), Base.length.(categorized_paths))
    return fig
end

function print_compactness_table1(depth::Integer, dim::Integer)
    padwidth = 10
    print(lpad("Length", padwidth)*lpad("Number", padwidth))
    for i in 0:9
        print(lpad("t=$i", padwidth))
    end
    print("\n")
    SAPs = topological_SAPs(depth, dim)
    for len in 3:depth
        categorized_SAPs = categorize_by_compactness(SAPs[len])
        print(lpad(len, padwidth)*lpad(length(SAPs[len]), padwidth))
        for category in categorized_SAPs
            print(lpad(length(category), padwidth))
        end
        print("\n")
    end
    print("\n")
end

function print_histogram(len::Integer, dim::Integer)
    SAPs = topological_SAPs(len, dim)
    categorized_SAPs = categorize_by_compactness(SAPs[len])
    plot_compactness_categorized(categorized_SAPs)
end

print_compactness_table1(16, 2)
print_histogram(16, 2)


