using Plots
using LaTeXStrings

h_range = LinRange(0.0, 1.5, 40)
Δτ_range = LinRange(0.001, 0.02, 20)

heatmap(h_range, Δτ_range, [(-log(tanh(Δτ * h)) / 2) / atanh(exp(- 2 * Δτ * h)) for Δτ in Δτ_range, h in h_range],
    xlabel = L"h", ylabel = L"\Delta \tau")

##

Δτ = 0.01
h_range = LinRange(0.0, 500, 40)

p = plot(h_range, atanh.(exp.(- 2 * Δτ * h_range)))
plot!(p, h_range, -log.(tanh.(Δτ * h_range)) / 2)

##

Δτ_range = LinRange(0.0, 1.0, 50)
h = 1

p = plot(Δτ_range, atanh.(exp.(- 2 * Δτ_range * h)))
plot!(p, Δτ_range, -log.(tanh.(Δτ_range * h)) / 2)