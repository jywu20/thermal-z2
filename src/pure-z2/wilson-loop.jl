using Random
using Statistics
using ProgressMeter
using Plots
using LaTeXStrings

include("utils.jl")
include("anisotropic_ising.jl")

##

function corr(ising_field)
    chain = ising_field[1, :]
    res = zeros(length(chain))
    for L in 1 : length(chain)
        res[L] = chain[1] * chain[L]
    end
    res
end

function wilson_loop(corr_funcs, corr_std)
    len = length(corr_funcs)
    perimeter = Array{Float64}(undef, (len - 1))
    area = Array{Float64}(undef, (len - 1))
    loop = Array{Float64}(undef, (len - 1))
    std_err = Array{Float64}(undef, (len - 1))
    for i in 2:len
        area[i - 1] = (i - 1) * (i - 1)
        perimeter[i - 1] = 4 * (i - 1)
        loop[i - 1] = abs(corr_funcs[i]^(i-1))
        std_err[i - 1] = (i - 1) * corr_funcs[i]^(i - 2) * corr_std[i]
    end
    (area=area, perimeter=perimeter, loop=loop, std_err=std_err)
end

function wilson_loop(; J = 1.0, h = 1.0, T = 1.0, n_sites = 10, n_steps = 100, n_heat = 10, n_sweep = 400)
    β = 1.0 / T
    dtau = β / n_steps
    j_tau = -log(tanh(dtau * h)) / 2

    kx = j_tau
    ky = J * dtau
    nx = n_steps
    ny = n_sites

    sim_params_heat = AnisotropicIsing2DSimParams(nx, ny, n_heat)
    ham_params = AnisotropicIsing2DHamParams(kx, ky)
    ising_field = init(sim_params_heat)
    run!(ham_params, sim_params_heat, ising_field)

    sim_params = AnisotropicIsing2DSimParams(nx, ny, n_sweep)
    
    corr_funcs, corr_err = binning(() -> run!(ham_params, sim_params, ising_field; observe = corr), n_bin)

    wilson_loop(corr_funcs, corr_err)
end

##
J = 1.0
h = 1.7
T = 0.01

n_steps = 100
n_sites = 20
n_heat = 500
n_bin = 10
n_sweep = 2000
dh = 0.05

area, perimeter, loop, std_err = wilson_loop(J = J, h = h, T = T, n_sites = n_sites, n_steps = n_steps, n_heat = n_heat, n_sweep = n_sweep)

cuthalf(arr) = arr[1 : Int64(n_sites / 2)]

area, perimeter, loop, std_err = map(cuthalf, (area, perimeter, loop, std_err))

##
scatter(area, log.(loop), yerror = std_err ./ loop, legend = false)

##
scatter(perimeter, log.(loop), yerror = std_err ./ loop, legend = false)
