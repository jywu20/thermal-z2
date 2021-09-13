using Plots
using Statistics
using ProgressMeter
using LaTeXStrings

include("ising-3d.jl")

function magnetization(;Jx, Jy, Jz, β, nx, ny, nz, n_bin, n_heat, n_sweep)
    sim_params_heat = AnisotropicIsingSimParams3DWolff(nx, ny, nz, n_heat)
    ham_params = AnisotropicIsingModel3D(β * Jx, β * Jy, β * Jz)
    ising_field = init_field_3d(Int64, sim_params_heat)
    run!(ising_field, ham_params, sim_params_heat)

    sim_params = AnisotropicIsingSimParams3DWolff(nx, ny, nz, n_sweep)
    mag, _ = binning(() -> run!(ising_field, ham_params, sim_params; observe = abs ∘ mean), n_bin)

    mag
end

##

T_range = LinRange(0.01, 6, 60)
J_range = LinRange(0.01, 4, 50)

progress = Progress(length(T_range) * length(J_range))
mag = zeros(length(T_range), length(J_range))

n_side = 10
n_bin = 10
n_heat = 100
n_sweep = 100

for (i, T) in enumerate(T_range)
    for (j, J) in enumerate(J_range)
        mag[i, j] = magnetization(
            Jx = J, Jy = J, Jz = J, β = 1 / T, 
            nx = n_side, ny = n_side, nz = n_side, 
            n_bin = n_bin, n_heat = n_heat, n_sweep = n_sweep)
        next!(progress)
    end
end

##

p = heatmap(J_range, T_range, mag, xlabel = L"J", ylabel = L"T")
J_range_analytical = LinRange(J_range[1], J_range[end], 50)
plot!(p, J_range_analytical, J_range_analytical / 0.2216544, legend = false)
ylims!(p, T_range[1], T_range[end])
xlims!(p, J_range[1], J_range[end])

##

J_one_index = argmin(@. abs(J_range - 1))

plot(T_range, mag[:, J_one_index])