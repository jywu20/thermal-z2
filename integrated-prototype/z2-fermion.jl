#region DQMC of fermions coupled to a Z2 gauge field.

# Note that since the fermions in the hopping Hamiltonian H_{hopping} = - t \sum_{⟨i, j⟩} σᶻ_{ij} c†_i c_j
# can be directly integrated out, the contribution of fermions in the discrete path integral is just the 
# determinant of B-matrices, and the "auxiliary field" here is just the Z2 gauge field itself.

"""
`A` is the type of the auxiliary field itself.

`M` is the type of matrices used in DQMC.
"""
abstract type AbstractFermionAuxField{L <: AbstractLattice, 
    A <: AbstractDiscretePathIntegralConfiguration, 
    F <:AbstractFloat} <: AbstractDiscretePathIntegralConfiguration end

"""
Convention: B[site1, site2, imaginary time]
"""
struct Z2SpinlessFermionAuxField{L <: AbstractLatticeWithPlaquattes, V, F} <: AbstractFermionAuxField{L, Z2GaugeFieldDPI{L, V}, F}
    σ::Z2GaugeFieldDPI{L, V}
    B::Array{F, 3}
end

struct Z2SpinlessFermionDQMC{F <: AbstractFloat, M <: AbstractMatrix}
    t::F
    Δτ::F
    β::F
    # TODO: Delta matrices
end

"""
B-matrices by definition, without any numerical acceleration.
"""
function B_def(model::Z2SpinlessFermionDQMC, σ::Z2GaugeFieldDPI, τ)
    Δτ = model.Δτ
    t = model.t
    n_site = site_number(σ.lattice)
    T = zeros(n_site, n_site)
    for b in bonds(lattice)
        i, j = bond_to_sites(lattice, b)
        T[i, j] = T[j, i] = σ[b, τ]
    end
    exp(Δτ * t * T)
end

function T_def(model::Z2SpinlessFermionDQMC, σ::Z2GaugeFieldDPI, τ)
    Δτ = model.Δτ
    t = model.t
    n_site = site_number(σ.lattice)
    T = zeros(n_site, n_site)
    for b in bonds(lattice)
        i, j = bond_to_sites(lattice, b)
        T[i, j] = T[j, i] = σ[b, τ]
    end
    - t * T
end

function B_at_bond(model::Z2SpinlessFermionDQMC{F, M}, σ::Z2GaugeFieldDPI, b, τ) where {F <: AbstractFloat, M <: AbstractMatrix}
    lattice = σ.lattice
    n_site = site_number(lattice)
    i, j = bond_to_sites(lattice, b)
    t = model.t
    Δτ = model.Δτ

    T = zeros(F, n_site, n_site)
    T[i, j] = T[j, i] = Δτ * t * σ[b, τ]
    M(exp(T))
end

B_at_bond(model::Z2SpinlessFermionDQMC, config::Z2SpinlessFermionAuxField, b, τ) = B_at_bond(model, config.σ, b, τ)
