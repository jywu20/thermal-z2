using Plots
using LinearAlgebra

include("../../z2-direct-prototype/lib.jl")

##
# Check whether the definitions in z2-direct-prototype\dqmc.jl is correct.

t = 1.0
β = 2.0
Δτ = β / n_steps
n_steps = 100
n_side = 4
n_sites = n_side^2
lattice = PeriodicSquareLattice2DWithBonds(n_side)
params = DeconfinedFermionsDQMCParams(t, n_side, n_steps, β)
σ = all_ones_discrete_path_integral_z2_gauge_field(lattice, n_steps)
τ = 1
T_real = z2_to_tight_binding_kinetic_matrix(σ, τ)
heatmap(T_real)

##
T_from_hopping = zeros(n_sites, n_sites)
for b in bonds(lattice) 
    T_from_hopping += bond_hopping_mat(σ, b, τ)
end
heatmap(T_from_hopping)

##
@show relative_err(T_from_hopping, T_real)

# We can say that the sum of the hopping terms are the correct kinetic Hamiltonian of the deconfined fermions.

##
# Define the B-matrices.

function B_τ(σ::DiscretePathIntegralZ2GaugeFieldPeriodicSquare2D, params::DeconfinedFermionsDQMCParams, τ::Int)
    Δτ = params.Δτ
    t = params.t
    exp(Δτ * t * z2_to_tight_binding_kinetic_matrix(σ, τ))
end

function B_β_τ(σ::DiscretePathIntegralZ2GaugeFieldPeriodicSquare2D, params::DeconfinedFermionsDQMCParams, τ::Int)
    n_τ = temporal_steps_count(σ)
    B = I
    for τ′ in τ + 1 : n_τ
        B = B_τ(σ, params, τ′) * B
    end
    B
end

function B_τ_0(σ::DiscretePathIntegralZ2GaugeFieldPeriodicSquare2D, params::DeconfinedFermionsDQMCParams, τ::Int)
    B = I
    for τ′ in 1 : τ
        B = B_τ(σ, params, τ′) * B
    end
    B
end

## 
# Verify whether these definitions agree with each other.

B_3_0_forward_one_step_all_ones = B_τ(σ, params, 3) * B_τ_0(σ, params, 2)
B_3_0_forward_all_ones = B_τ_0(σ, params, 3)
@show relative_err(B_3_0_forward_all_ones, B_3_0_forward_one_step_all_ones)

#
σ[1, 2, 1] *= -1
σ[2, 6, 2] *= -1
B_5_0_forward_one_step_flipped_1 = B_τ(σ, params, 5) * B_τ_0(σ, params, 4)
B_5_0_forward_flipped_1 = B_τ_0(σ, params, 5)
@show relative_err(B_5_0_forward_one_step_flipped_1, B_5_0_forward_one_step_flipped_1)

##
σ[4, 16, 5] *= -1
σ[3, 2, 4] *= -1
B_5_0_forward_one_step_flipped_2 = B_τ(σ, params, 5) * B_τ_0(σ, params, 4)
B_5_0_forward_flipped_2 = B_τ_0(σ, params, 5)
@show relative_err(B_5_0_forward_one_step_flipped_2, B_5_0_forward_one_step_flipped_2)

##
σ[4, 16, 5] *= -1
σ[3, 2, 4] *= -1
σ[3, 4, 5] *= -1
B_10_5_backward_one_step_flipped_2 = B_β_τ(σ, params, 5) * B_τ(σ, params, 5)
B_10_5_flipped_2 = B_β_τ(σ, params, 4)
@show relative_err(B_10_5_backward_one_step_flipped_2, B_10_5_flipped_2)

# ?? The error is not exactly zero. Maybe some problem, maybe not.
# It seems that there is no real problem. The small error seems to be caused by the fact that 5 is far from n_steps.

##
# Complete Trotter decomposition.

function B_τ_trotter(σ::DiscretePathIntegralZ2GaugeFieldPeriodicSquare2D, params::DeconfinedFermionsDQMCParams, τ::Int)
    lattice = σ.lattice 
    t = params.t
    Δτ = params.Δτ
    B = I
    for b in bonds(lattice)
        B = exp(Δτ * t * bond_hopping_mat(σ, b, τ)) * B
    end
    B
end

@show relative_err(B_τ_trotter(σ, params, 10), B_τ(σ, params, 10))