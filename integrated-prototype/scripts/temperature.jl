using Plots

##
Δτ = 0.01

T_range = LinRange(0.5, 5, 50)
β_range = 1 ./ T_range
n_τ_range = Int.(round.(β_range / Δτ))

β′_range = Δτ * n_τ_range
T′_range = 1 ./ β′_range

plot(T_range, T′_range, legend = false, xlabel = "input T", ylabel = "actual T")