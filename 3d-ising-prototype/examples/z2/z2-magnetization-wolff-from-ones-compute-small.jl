using Statistics
using JLD

include("../../src/ising-3d.jl")

println("Package initialization completed.")

function magnetization(;J, h, β, n_side, n_steps, n_bin, n_heat, n_sweep)
    sim_params_heat = TransverseFieldIsingSimParams2DWolff(β, n_steps, n_side, n_heat)
    # J_ising = h_z2, h_ising = J_z2
    ham_params = TransverseFieldIsingModel2D(h, J)
    ising_field = init_field_3d_ones(Int64, AnisotropicIsingSimParams3DWolff(n_side, n_side, n_steps, n_sweep))
    run!(ising_field, ham_params, sim_params_heat)

    sim_params = TransverseFieldIsingSimParams2DWolff(β, n_steps, n_side, n_sweep)
    mag, _ = binning(() -> run!(ising_field, ham_params, sim_params; 
        observe = abs ∘ mean ∘ (x -> x[:, :, 1])), n_bin)

    mag
end

T_steps = 80
T_range = LinRange(0.1, 4.1, T_steps)
h_steps = 80
h_range = LinRange(0.1, 3.1, h_steps)

J = 1.0

mag = zeros(length(T_range), length(h_range))

n_side = 10
n_steps = 20
n_bin = 10
n_heat = 2000
n_sweep = 1000

output_text_name = "output"
output_jld_name = "data.jld"

open(output_text_name, "w") do file
    println(file, "n_side   =   $n_side")
    println(file, "n_steps  =   $n_steps")
    println(file, "n_bin    =   $n_bin")
    println(file, "n_heat   =   $n_heat")
    println(file, "n_sweep  =   $n_sweep")
    println(file, "T_min    = $(T_range[1]),    T_max = $(T_range[end]),   steps = $T_steps")
    println(file, "h_min    = $(h_range[1]),    T_max = $(h_range[end]),   steps = $h_steps")
    println(file, "J        =   $J")
end

println("Monte Carlo calculation starts.")

for (i, T) in enumerate(T_range)
    for (j, h) in enumerate(h_range)
        mag[i, j] = magnetization(
            J = J, h = h, 
            β = 1 / T, n_side = n_side, n_steps = n_steps,
            n_bin = n_bin, n_heat = n_heat, n_sweep = n_sweep
        )

        open(output_text_name, "a") do file
            println(file, "T = $T, h = $h, mag = $(mag[i, j])")
        end

        println("T = $T, h = $h completed.")
    end
end
