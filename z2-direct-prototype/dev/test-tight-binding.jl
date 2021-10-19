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

for τ in [1, 1, 1, 1, 1]
    T_origin = z2_to_tight_binding_kinetic_matrix(σ, τ)
    U_flipped_exact = exp(- Δτ * T_origin)
    U_flipped_estimated = exp(- Δτ * T_origin) 
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

##

τ = 2
T_real = z2_to_tight_binding_kinetic_matrix(σ, τ)
T_from_plaquattes = sum([plaquatte_hamiltonian(σ, p, τ) for p in sites(lattice)]) / 2

relative_err(T_real, T_from_plaquattes)

##
T_A, T_B = checkboard_decomposition_hamiltonians(lattice, σ, τ)

(T_A + T_B) == T_real

##

Δτ = 0.001
relative_err(exp(- Δτ * T_real), exp(- Δτ * T_A) * exp(- Δτ * T_B))

##

Δτ = 0.001
errors = Float64[]

for τ in 1 : n_τ
    for i in sites(lattice)
        for j in nearest_neighbors(lattice, i)
            if rand() < 0.4
                σ[i, j, τ] *= -1
                T_real = z2_to_tight_binding_kinetic_matrix(σ, τ)
                T_A, T_B = checkboard_decomposition_hamiltonians(lattice, σ, τ)
            end
            push!(errors, relative_err(exp(- Δτ * T_real), exp(- Δτ * T_A) * exp(- Δτ * T_B)))
        end
    end
end

plot(errors)

##

Δτ = 0.01
errors = Float64[]

n_sites = length(sites(lattice))

sub_A, sub_B = checkboard_decomposition_bonds(lattice)

for τ in 1 : 1
    for b in [1]#bonds(lattice)
        T_real = z2_to_tight_binding_kinetic_matrix(σ, τ)
        T_A, T_B = checkboard_decomposition_hamiltonians(lattice, σ, τ)
        B_A = exp(- Δτ * T_A)
        B_B = exp(- Δτ * T_B)
        B_real = exp(- Δτ * T_real)
        if rand() < 1.0
            σ[b, τ] *= -1
            i, j = bond_to_sites(lattice, b)
            Δ = exp(- 2Δτ * Tij(n_sites, i, j))
            if b in sub_A
                B_A = Δ * B_A
            else
                B_B = Δ * B_B
            end
        end
        push!(errors, relative_err(B_real, B_A * B_B))
    end
end

#plot(errors, legend = false)

##
# Problem: in the following definition, Δ does not commute with hopping terms in the plaquatte of which one bond is flipped.
# Maybe we should use a 4 × 4 Δ to lessen the error.
τ = 1
b = 1
T_real = z2_to_tight_binding_kinetic_matrix(σ, τ)
T_A, T_B = checkboard_decomposition_hamiltonians(lattice, σ, τ)
B_A = exp(- Δτ * T_A)
B_B = exp(- Δτ * T_B)
B_real = exp(- Δτ * T_real)
σ[b, τ] *= -1
i, j = bond_to_sites(lattice, b)
Δ = exp(- 2Δτ * Tij(n_sites, i, j))
if b in sub_A
    B_A = Δ * B_A
else
    B_B = Δ * B_B
end

@show relative_err(B_real, B_A * B_B)

##

plaquatte_shared_bond(lattice, 1)

##

p = 1
τ = 1
p_bonds = plaquatte(lattice, p)
@show p_bonds
heatmap(plaquatte_hamiltonian(σ, p, τ))
@show map(x -> σ[x, τ], p_bonds)

##

τ = 1
b = 1
T_A, T_B = checkboard_decomposition_hamiltonians(lattice, σ, τ)
B_A = exp(- Δτ * T_A)
B_B = exp(- Δτ * T_B)
σ[b, τ] *= -1
i, j = bond_to_sites(lattice, b)
Δ = exp(- 2Δτ * Tij(n_sites, i, j) * σ[b, τ])
if b in sub_A
    B_A = Δ * B_A
else
    B_B = Δ * B_B
end
T_real = z2_to_tight_binding_kinetic_matrix(σ, τ)
B_real = exp(- Δτ * T_real)

@show relative_err(B_real, B_A * B_B)