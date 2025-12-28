using SantaLucia
using CSV
using DataFrames
using GLMakie

data_filepath = "data/Watson-Crick-NN-Parameters.csv"
data_inner = CSV.read(data_filepath, DataFrame)

println("-----------------------------------------------------------")

#CGTTGA

sequence_string = "CGTTGA"
C = 1 

sequence = Sequence(sequence_string) 

potentials = thermodynamic_potentials(sequence, data_inner)

H = potentials[1]
S = potentials[2]

Tm = melting_temperature(potentials...,C)

print("Strand = ", sequence_string, "\n")

print("ΔH = ")
print(round(H/1000, digits = 1))
print(" kcal/mol \n")

print("ΔS = ")
print(round(S, digits = 1))
print(" cal/mol \n")

print("Tm = ")
print(Tm-273.15)
print(" C° (")
print(Tm)
print(" K)\n")

println("-----------------------------------------------------------")

Ts = Vector(0:1:600)
Xs = Vector{Float64}(undef, length(Ts))

for i in 1:length(Xs)
    Xs[i] = melting_curve(1/Ts[i], H, S)
end

fig = Figure(size=(800, 600), fontsize = 25)
ax = Axis(
    fig[1,1],
    xlabel = "Temeperature (K)",
    ylabel = "X"
    )
xlims!(ax, 0, 600)
lines!(ax, Ts, Xs, linewidth = 3)
display(fig)

