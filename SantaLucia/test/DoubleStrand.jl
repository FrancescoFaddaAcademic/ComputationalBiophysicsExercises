using SantaLucia
using CSV
using DataFrames
using GLMakie

data_filepath = "data/Watson-Crick-NN-Parameters.csv"
inner_data = CSV.read(data_filepath, DataFrame)

println("-----------------------------------------------------------")

#CGTTGA

sequence_string1 = "CGTTGA"
sequence_string2 = "TCAACG"

sequence1 = Sequence(sequence_string1)
sequence2 = Sequence(sequence_string2)

if !are_complementary(sequence1, sequence2)
    error("The two sequences are not complementary")
end

C = 1 

potentials = sequence_potentials(sequence1, inner_data)

H = potentials[1]
S = potentials[2]

Tm = melting_temperature(potentials...,C)

print("Strand = ", sequence_string1, "\n")

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

Ts = Vector(0:0.0001:600)
Xs = Vector{Float64}(undef, length(Ts))
δ = 0.001

for i in 1:length(Xs)
    Xs[i] = melting_curve(1/Ts[i], H, S)
end

particular = findall(x -> (1-δ)>x>(0+δ), Xs)
fig = Figure(size=(1280, 1090), fontsize = 20)

ax = Axis(
    fig[1,1],
    xlabel = "Temeperature (K)",
    ylabel = "X",
    yticks = 0:0.1:1,
    xticks = WilkinsonTicks(10, k_min=9),
    )
ylims!(ax, -0.1, 1.1)

scatter!(ax, Ts[particular], Xs[particular], markersize = 5)
display(fig)
