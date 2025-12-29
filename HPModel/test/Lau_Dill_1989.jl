using HPModel
using GLMakie
using LaTeXStrings

#This test reproduces Figure8, Figure12 and Figure13 from the article [Lau-Dill 1989] 

function plot_g(g::Vector{Int}, ttl)
    bins = Base.length(g)
    binwidth = 1
    fig = Figure(size=(800, 600), fontsize = 15)
    ax1 = Axis(
        fig[1, 1],
        title=ttl,
        xlabel = L"g(s)", 
        ylabel = "N",
        xgridvisible=false,
        ygridvisible=true,
        xticklabelsize = 20,
        yticklabelsize = 20,
        xticks = round.(0:ceil(bins/100)*10:bins, digits = 2),
    )
    ax2 = Axis(
        fig[2, 1],
        xlabel = L"g(s)", 
        ylabel = "N",
        xgridvisible=false,
        ygridvisible=true,
        xticklabelsize = 20,
        yticklabelsize = 20,
        xticks = round.(0:1:20, digits = 2),
    ) 
    xlims!(ax1, 1-binwidth/2, bins+binwidth/2)
    ylims!(ax1, 0, nothing)
    barplot!(ax1, g)
    xlims!(ax2, 1-binwidth/2, 20+binwidth/2)
    ylims!(ax2, 0, nothing)
    barplot!(ax2, g)
    display(fig)
    #return fig
end

dim = 2

len = 10
title = "Figure 8: Distribution of the native-state properties for the full sequence space (n = 10)"
SAPs = unique_SAPs(len, dim)
sequences = all_binary_sequences(len)
self_adjacency_triangles = path_self_adjacency_triangle.(SAPs)
g = calculate_g(sequences, self_adjacency_triangles, multiplicative_binary_interaction_model)
plot_g(g, title)

len = 10
title = "Figure 12: Distribution of the native-state properties for the most compact configurations (n = 10)"
SAPs = unique_SAPs(len, dim)
categorized_SAPs = categorize_by_compactness(SAPs)
self_adjacency_triangles = path_self_adjacency_triangle.(categorized_SAPs[length(categorized_SAPs)])
sequences = all_binary_sequences(len)
g = calculate_g(sequences, self_adjacency_triangles, multiplicative_binary_interaction_model)
plot_g(g, title)

len = 13
title = "Figure 13: Distribution of the native-state properties \n for the most compact configurations of 200 random sequences (n = 13)"
SAPs = unique_SAPs(len, dim)
categorized_SAPs = categorize_by_compactness(SAPs)
self_adjacency_triangles = path_self_adjacency_triangle.(categorized_SAPs[length(categorized_SAPs)])
sequences = random_binary_sequences(len, 200)
g = calculate_g(sequences, self_adjacency_triangles, multiplicative_binary_interaction_model)
plot_g(g, title)


