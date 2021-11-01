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
    U_L::Array{F, 3}
    U_R::Array{F, 3}
    D_L::Array{F, 3}
    D_R::Array{F, 3}
    V_L::Array{F, 3}
    V_R::Array{F, 3}
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

function B_τ_0_def(model::Z2SpinlessFermionDQMC, σ::Z2GaugeFieldDPI, τ)
    n_site = site_number(σ.lattice)
    B = Matrix(I, n_site, n_site)
    for τ′ in 1 : τ
        B = B_def(model, σ, τ′) * B        
    end
    B
end

function B_β_τ_def(model::Z2SpinlessFermionDQMC, σ::Z2GaugeFieldDPI, τ)
    n_site = site_number(σ.lattice)
    n_τ = σ.n_τ
    B = Matrix(I, n_site, n_site)
    for τ′ in n_τ : τ + 1
        B = B * B_def(model, σ, τ′)
    end
    B
end

function G_τ_def(model::Z2SpinlessFermionDQMC, σ::Z2GaugeFieldDPI, τ)
    inv(I + B_τ_0_def(model, σ, τ) * B_β_τ_def(model, σ, τ))
end

function T_def(model::Z2SpinlessFermionDQMC, σ::Z2GaugeFieldDPI, τ)
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

function B_prod_bonds(model::Z2SpinlessFermionDQMC{F, M}, σ::Z2GaugeFieldDPI, τ) where {F <: AbstractFloat, M <: AbstractMatrix}
    lattice = σ.lattice
    map(b -> B_at_bond(model, σ, b, τ), bonds(lattice)) |> prod
end

function Δ_mat(model::Z2SpinlessFermionDQMC{F, M}, σ::Z2GaugeFieldDPI, b, τ) where {F <: AbstractFloat, M <: AbstractMatrix}
    Δτ = model.Δτ
    t = model.t
    lattice = σ.lattice
    n_site = site_number(lattice)
    i, j = bond_to_sites(lattice, b)
    Δ = zeros(F, n_site, n_site)
    Δ[i, i] = Δ[j, j] = cosh(2 * Δτ * t * σ[b, τ]) - 1
    Δ[i, j] = Δ[j, i] = - sinh(2 * Δτ * t * σ[b, τ])
    Δ
end