abstract type AbstractModel{F <: Float64} end

local_hamiltonian_cluster(model::M) where{M <: AbstractModel{F}, F <: Float64}
    = @error "Model $M"