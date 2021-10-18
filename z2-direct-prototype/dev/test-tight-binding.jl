using Plots
using Statistics
using LinearAlgebra

include("../lib.jl")

# Tiny example

T0 = [
    0 1 1 0 ;
    1 0 0 1 ;
    1 0 0 1 ; 
    0 1 1 0
]

# Change the hopping term between site 1 and 2
Tp = [
    0 -1 1 0 ;
    -1 0 0 1 ;
    1 0 0 1 ; 
    0 1 1 0
]

# Calculate Δ
Δτ = 0.001
Δ = exp(- Δτ * Tp) * exp(Δτ * T0) - I

Δ_approx_1 = exp(- Δτ * (Tp - T0)) - I

relative_err(Δ, Δ_approx_1)

##

lattice = PeriodicSquareLattice2DWithBonds(4)
T = tight_binding_kinetic_matrix_periodic(lattice)

Tp = copy(T)
Tp[1, 2] *= -1
Tp[2, 1] *= -1

Δτ = 0.001
Δ = exp(- Δτ * Tp) * exp(Δτ * T) - I

Δ_approx_1 = exp(- Δτ * (Tp - T)) - I

relative_err((I + Δ) * exp(- Δτ * T), (I + Δ_approx_1) * exp(- Δτ * T))

##
n_τ = 20
σ = all_ones_discrete_path_integral_z2_gauge_field(lattice, n_τ)
T_origin = z2_to_tight_binding_kinetic_matrix(σ, 1)

σ[1, 5, 1] *= -1
T_change_1 = z2_to_tight_binding_kinetic_matrix(σ, 1)
Δ_1 = exp(- Δτ * (T_change_1 - T_origin))

σ[4, 8, 1] *= -1
T_change_2 = z2_to_tight_binding_kinetic_matrix(σ, 1)
Δ_2 = exp(- Δτ * (T_change_2 - T_change_1))

σ[7, 11, 1] *= -1
T_change_3 = z2_to_tight_binding_kinetic_matrix(σ, 1)
Δ_3 = exp(- Δτ * (T_change_3 - T_change_2))

relative_err(Δ_3 * Δ_2 * Δ_1 * exp(- Δτ * T_origin), exp(- Δτ * T_change_3))

##

Δτ = 0.001

τ = 1
n_τ = 20
σ = all_ones_discrete_path_integral_z2_gauge_field(lattice, n_τ)
T_origin = z2_to_tight_binding_kinetic_matrix(σ, τ)
U_flipped_exact = exp(- Δτ * T_origin)
U_flipped_estimated = exp(- Δτ * T_origin) 

function Δ_mat(n, i, j)
    Δ = zeros(n, n)
    Δ[i, j] = 1
    Δ[j, i] = 1
    Δ    
end

function Tij(n_sites, i, j)
    res = zeros(n_sites, n_sites)
    res[i, j] = 1
    res[j, i] = 1
    res
end

sites_length = length(sites(lattice))
errors = Float64[]

for τ in 1 : n_τ
    for i in sites(lattice)
        for j in nearest_neighbors(lattice, i)
            if rand() < 0.5
                σ[i, j, τ] *= -1
                U_flipped_estimated = exp(- 2Δτ * Δ_mat(sites_length, i, j) * σ[i, j, τ]) * U_flipped_estimated
                U_flipped_exact = exp(- Δτ * z2_to_tight_binding_kinetic_matrix(σ, τ))
            end
            push!(errors, relative_err(U_flipped_exact, U_flipped_estimated))
        end
    end
end

plot(errors, legend = false) 

##
checkboard_decomposition(lattice)