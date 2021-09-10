using Random
using Statistics
using ProgressMeter
using Plots
using LaTeXStrings

include("utils.jl")
include("anisotropic_ising.jl")

##

# We don't use `observe.jl` because it averages over all spins in the 2D Ising model, but since one dimension
# in that 2D Ising model is actually the time dimension, we should only average over spins where τ = β
function magnetization(;kx, ky, nx, ny, n_bin, n_heat, n_sweep)
    sim_params_heat = AnisotropicIsing2DSimParams(nx, ny, n_heat)
    ham_params = AnisotropicIsing2DHamParams(kx, ky)
    ising_field = init(sim_params_heat)
    run!(ham_params, sim_params_heat, ising_field)

    sim_params = AnisotropicIsing2DSimParams(nx, ny, n_sweep)
    mag, _ = binning(() -> run!(ham_params, sim_params, ising_field; observe = abs ∘ mean ∘ (x -> x[1, :])), n_bin)

    mag
end

##
J = 1.0
h = 0.0

n_steps = 100
n_sites = 10
n_heat = 10
n_bin = 10
n_sweep = 400
dh = 0.05

T_range = LinRange(0.05, 4, 20)
h_range = LinRange(0.01, 2.5, 40)
mag_to_T_h = Array{Float64}(undef, length(T_range), length(h_range))
mag_to_T_h_dh = Array{Float64}(undef, length(T_range), length(h_range))
progress = Progress(length(T_range) * length(h_range))

for (i, T) in enumerate(T_range), (j, h) in enumerate(h_range)
    β = 1.0 / T
    dtau = β / n_steps
    j_tau = -log(tanh(dtau * h)) / 2
    mag_to_T_h[i, j] = magnetization(kx=j_tau, ky=J * dtau, nx=n_steps, ny=n_sites, n_bin=n_bin, n_heat=n_heat, n_sweep=n_sweep)

    j_tau_dh = -log(tanh(dtau * (h + dh))) / 2
    mag_to_T_h_dh[i, j] = magnetization(kx=j_tau_dh, ky=J * dtau, nx=n_steps, ny=n_sites, n_bin=n_bin, n_heat=n_heat, n_sweep=n_sweep)

    next!(progress)
end

##

heatmap(h_range, T_range, mag_to_T_h, 
    xlabel = L"h", 
    ylabel = L"T",
    title = "Magnetization"
)

##

heatmap(h_range, T_range, (mag_to_T_h_dh - mag_to_T_h) / dh,
    xlabel = L"h", 
    ylabel = L"T",
    title = "Magnetic susceptibility"
)
