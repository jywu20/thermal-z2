#region field configuration for path integral of transverse field Ising model

"""
`L` is the lattice type and `D` is the data type of Ising degrees of freedom.
"""
struct DiscretePathIntegralIsingField{L <: AbstractLattice, D}
    n_steps::Int
    lattice::L
    data::Array{D, 2}
end

getindex(ising_field::DiscretePathIntegralIsingField, idx...) = ising_field.data[idx...]
getindex(ising_field::DiscretePathIntegralIsingField, i::Int, τ::Int) = ising_field.data[i, τ]
setindex!(ising_field::DiscretePathIntegralIsingField, value, idx...) = ising_field.data[idx...] = value
setindex!(ising_field::DiscretePathIntegralIsingField, value, i::Int, τ::Int) = ising_field.data[i, τ] = value

function all_ones_discrete_path_integral_ising_field(
    ::Type{D}, lattice::L, n_steps::Int
) where {L <: AbstractLattice}
    n_sites = length(sites(lattice))
    data = ones(D, n_sites, n_steps)
    DiscretePathIntegralIsingField(n_steps, lattice, data)
end 

all_ones_discrete_path_integral_ising_field(lattice::L, n_steps::Int) where {L <: AbstractLattice} = 
    all_ones_discrete_path_integral_ising_field(Int, lattice, n_steps)

const DiscretePathIntegralIsingFieldPerdiodicSquare2D{D} = 
    DiscretePathIntegralIsingField{PeriodicSquareLattice2DWithBonds, D}

#endregion

#region the TFIM model: Metropolis algorithm

