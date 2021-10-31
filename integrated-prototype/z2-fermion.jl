#region DQMC of fermions coupled to a Z2 gauge field.

# Note that since the fermions in the hopping Hamiltonian H_{hopping} = - t \sum_{⟨i, j⟩} σᶻ_{ij} c†_i c_j
# can be directly integrated out, the contribution of fermions in the discrete path integral is just the 
# determinant of B-matrices, and the "auxiliary field" here is just the Z2 gauge field itself.

struct Z2SpinlessFermionDQMC{L <: AbstractLatticeWithPlaquattes{Int, Int, Int}, V, F <:AbstractFloat} <: AbstractDiscretePathIntegralConfiguration
    σ::Z2GaugeFieldDPI{L, V}
    B::Array{F} # TODO:size
end


