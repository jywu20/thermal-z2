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
    M <: AbstractMatrix,
    F <:AbstractFloat} <: AbstractDiscretePathIntegralConfiguration end

"""
Convention: B[site1, site2, , imaginary time]
"""
struct Z2SpinlessFermionAuxField{L <: AbstractLatticeWithPlaquattes, V, M, F} <: AbstractFermionAuxField{L, Z2GaugeFieldDPI{L, V}, M, F}
    σ::Z2GaugeFieldDPI{L, V}
    B::Array{F, 3}
end

struct Z2SpinlessFermionDQMC{F <: AbstractFloat}
    t::F
    Δτ::F
    β::F
end

"""
B-matrices by definition, without any numerical acceleration.
"""
function B_def(model::Z2SpinlessFermionDQMC, σ::Z2GaugeFieldDPI, τ::Int)
    Δτ = model.Δτ
    t = model.t
    exp(Δτ * t * σ[:, :, τ])
end
