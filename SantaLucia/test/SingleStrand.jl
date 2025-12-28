using SantaLucia
using CSV
using DataFrames
using GLMakie

println("-----------------------------------------------------------")

sequence =  "TTTTCCCCTTTTCCCAAAACCCCAAAA"
structure = "((((....((((...))))....))))"

strand = SingleStrand(sequence, structure)

C = 1

hairPins = Vector{HairPin}([])

extract_hairpins!(strand, hairPins)

Δ = [0.0, 0.0]

hairpin_data = CSV.read("data/Hairpins-Turner-Parameters.csv", DataFrame)
internal_loop_data = CSV.read("data/Internal-Loops-Turner-Parameters.csv", DataFrame)
bulge_data = CSV.read("data/Bulges-Turner-Parameters.csv", DataFrame)
inner_data = CSV.read("data/Watson-Crick-NN-Parameters.csv", DataFrame)

for hairPin in hairPins
    global Δ += loop_thermodynamic_contribution(hairPin, hairpin_data, internal_loop_data, bulge_data)
    global Δ += sequence_potentials(hairPin.paired_sequence, inner_data)
end

H = Δ[1]
S = Δ[2]

Tm = melting_temperature(Δ...,C)

print("Sequence = ", sequence, "\n")
print("Structure = ", structure, "\n")

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


Ts = Vector(0:0.0001:360)
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