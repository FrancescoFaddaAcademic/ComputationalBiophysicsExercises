using HPModel
using LaTeXStrings
using GLMakie

#This test reproduces Table1 and Figure2 from the article [Chan-Dill 1989] 

function plot_compactedness_categorized(categorized_paths::Vector{Vector{Path}})
    bins = Base.length(categorized_paths)
    binwidth = 1/(bins-1)
    fig = Figure(size=(800, 600), fontsize = 25)
    ax = Axis(
        fig[1, 1], 
        #title = "Number of sequences by compactedness", 
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

function print_compactness_table(depth::Int)
    padwidth = 10
    print(lpad("Length", padwidth)*lpad("Number", padwidth))
    for i in 0:9
        print(lpad("t=$i", padwidth))
    end
    print("\n")
    dim = 2
    for len in 3:depth
        SAPs = unique_SAPs(len, dim)
        print(lpad(len, padwidth)*lpad(length(SAPs), padwidth))
        
        for category in categorized_SAPs
            print(lpad(length(category), padwidth))
        end
        print("\n")
    end
    print("\n")
end

function print_histogram()
    len = 16
    dim = 2
    SAPs = unique_SAPs(len, dim)
    categorized_SAPs = categorize_by_compactedness(SAPs)
    plot_compactedness_categorized(categorized_SAPs)
end

print_compactness_table(16)
print_histogram()

