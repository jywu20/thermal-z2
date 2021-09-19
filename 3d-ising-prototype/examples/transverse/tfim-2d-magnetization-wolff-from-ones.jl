using Plots
using Statistics
using ProgressMeter
using LaTeXStrings

include("../../src/ising-3d.jl")

function magnetization(;J, h, β, n_side, n_steps, n_bin, n_heat, n_sweep)
    sim_params_heat = TransverseFieldIsingSimParams2DWolff(β, n_steps, n_side, n_heat)
    ham_params = TransverseFieldIsingModel2D(J, h)
    ising_field = init_field_3d_ones(Int64, AnisotropicIsingSimParams3DWolff(n_side, n_side, n_steps, n_sweep))
    run!(ising_field, ham_params, sim_params_heat)

    sim_params = TransverseFieldIsingSimParams2DWolff(β, n_steps, n_side, n_sweep)
    mag, _ = binning(() -> run!(ising_field, ham_params, sim_params; 
        observe = abs ∘ mean ∘ (x -> x[:, :, 1])), n_bin)

    mag
end

##

T_steps = 60
T_range = LinRange(0.01, 6, T_steps)
ΔT = (T_range[end] - T_range[1]) / T_steps

h_steps = 50
h_range = LinRange(0.01, 4, h_steps)
Δh = (h_range[end] - h_range[1]) / h_steps

J = 1.0

progress = Progress(length(T_range) * length(h_range))
mag = zeros(length(T_range), length(h_range))

n_side = 10
n_steps = 20
n_bin = 10
n_heat = 100
n_sweep = 100

for (i, T) in enumerate(T_range)
    for (j, h) in enumerate(h_range)
        mag[i, j] = magnetization(
            J = J, h = h, 
            β = 1 / T, n_side = n_side, n_steps = n_steps,
            n_bin = n_bin, n_heat = n_heat, n_sweep = n_sweep
        )
        next!(progress)
    end
end

##

p = heatmap(h_range, T_range, mag, xlabel = L"h", ylabel = L"T")
ylims!(p, T_range[1], T_range[end])
xlims!(p, h_range[1], h_range[end])

##

mag_sus = [abs(mag[i, j + 1] - mag[i, j]) / Δh for i in 1 : length(T_range), j in 1 : length(h_range) - 1]
p = heatmap(h_range[1 : end - 1], T_range, mag_sus, xlabel = L"h", ylabel = L"T")
ylims!(p, T_range[1], T_range[end])
xlims!(p, h_range[1], h_range[end - 1])