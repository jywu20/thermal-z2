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
type parameters:
- `L`: the lattice.
- `V`: type of field value at each bond, for example `Int`
- `F`: the type of float used, for example `Float64` or `Float32`. 
  Note that `F` must support functions like `sin` or `exp`.
- `S`: a three-dimensional array type used to store auxiliary matrices, for example `Array{Float64, 3}`. 
The indexing convention: [site1, site2, imaginary time]
"""
struct Z2SpinlessFermionAuxField{
    L <: AbstractLatticeWithPlaquattes, V, F <: AbstractFloat, S <: AbstractArray
} <: AbstractFermionAuxField{L, Z2GaugeFieldDPI{L, V}, F}
    σ::Z2GaugeFieldDPI{L, V}
    U_L::S
    U_R::S
    D_L::S
    D_R::S
    V_L::S
    V_R::S
    G::S
end

struct Z2SpinlessFermionDQMC{F <: AbstractFloat, M <: AbstractMatrix}
    t::F
    Δτ::F
    β::F
    n_wrap::Int
    Δ_σ_1::Array{M, 1}
    Δ_σ_m1::Array{M, 1}
end

function Z2SpinlessFermionDQMC(::Type{F}, ::Type{M}, σ::Z2GaugeFieldDPI, n_wrap::Int, t::F, Δτ::F, β::F) where {
    F <: AbstractFloat, M <: AbstractMatrix}
    lattice = σ.lattice
    n_site = site_number(lattice)

    Δ_σ_1 = Array{M, 1}(undef, bond_number(lattice))
    Δ_σ_m1 = Array{M, 1}(undef, bond_number(lattice))
    for b in bonds(lattice)
        i, j = bond_to_sites(lattice, b)

        Δ = zeros(F, n_site, n_site)
        Δ[i, i] = Δ[j, j] = cosh(2 * Δτ * t) - 1
        Δ[i, j] = Δ[j, i] = - sinh(2 * Δτ * t)
        Δ_σ_1[b] = Δ

        Δ = zeros(F, n_site, n_site)
        Δ[i, i] = Δ[j, j] = cosh(- 2 * Δτ * t) - 1
        Δ[i, j] = Δ[j, i] = - sinh(- 2 * Δτ * t)
        Δ_σ_m1[b] = Δ
    end

    n_τ = time_step_number(σ)
    if n_τ % n_wrap != 0
        error("n_wrap must be a divisor of n_τ.")
    end

    Z2SpinlessFermionDQMC{F, M}(t, Δτ, β, n_wrap, Δ_σ_1, Δ_σ_m1)
end

function Δ_mat(model::Z2SpinlessFermionDQMC{F, M}, σ::Z2GaugeFieldDPI{L, V}, b, τ) where {
    L <: AbstractLattice, V, F <: AbstractFloat, M <: AbstractMatrix
}
    if σ[b, τ] == one(V) 
        model.Δ_σ_1[b]
    else
        model.Δ_σ_m1[b]
    end
end

include("defs.jl")

include("bond.jl")
