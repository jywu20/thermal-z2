import Base.getindex
import Base.setindex!
import Base.copy
import Base.+
import Base.-

#region field configuration for path integral of z2 gauge fields

"""
`L` is the lattice type and `D` is the data type of Z2 gauge degrees of freedom. 
Indices of `data` are bond indices of `lattice`, or site indices of `bond_lattice(lattice)`.
"""
struct DiscretePathIntegralZ2GaugeField{L <: AbstractLattice, D}
    n_steps::Int
    lattice::L
    data::Array{D, 2}
end

getindex(z2field::DiscretePathIntegralZ2GaugeField, idx...) = z2field.data[idx...]
getindex(z2field::DiscretePathIntegralZ2GaugeField, i::Int64, j::Int64, τ::Int64) = z2field.data[bond(z2field.lattice, i, j)..., τ]
setindex!(z2field::DiscretePathIntegralZ2GaugeField, value, idx...) = z2field.data[idx...] = value
setindex!(z2field::DiscretePathIntegralZ2GaugeField, value, i::Int64, j::Int64, τ::Int64) = z2field.data[bond(z2field.lattice, i, j)..., τ] = value

copy(z2field::DiscretePathIntegralZ2GaugeField) = DiscretePathIntegralZ2GaugeField(z2field.n_steps, z2field.lattice, copy(z2field.data))
# +(σ1::DiscretePathIntegralZ2GaugeField, σ2::DiscretePathIntegralZ2GaugeField) = σ1.data + σ2.data
# -(σ1::DiscretePathIntegralZ2GaugeField, σ2::DiscretePathIntegralZ2GaugeField) = σ1.data - σ2.data

function all_ones_discrete_path_integral_z2_gauge_field(
    ::Type{D}, lattice::L, n_steps::Int
) where {D, L <: AbstractLattice}
    n_sites = length(sites(bond_lattice(lattice)))
    data = ones(D, n_sites, n_steps)
    DiscretePathIntegralZ2GaugeField{L, D}(n_steps, lattice, data)
end

all_ones_discrete_path_integral_z2_gauge_field(lattice::L, n_steps::Int) where {L <: AbstractLattice} = 
    all_ones_discrete_path_integral_z2_gauge_field(Int, lattice, n_steps)

const DiscretePathIntegralZ2GaugeFieldPeriodicSquare2D{D} = 
    DiscretePathIntegralZ2GaugeField{PeriodicSquareLattice2DWithBonds, D}

# function DiscretePathIntegralZ2GaugeFieldPeriodicSquare2D(::Type{D}, lattice::PeriodicSquareLattice2DWithBonds, n_steps::Int) where {D} 
#     n_sites = length(sites(lattice))
#     DiscretePathIntegralZ2GaugeFieldPeriodicSquare2D(lattice, Array{D, 2}(undef, n_sites, n_steps))
# end

#endregion

#region the Z2 model: Metropolis algorithm

struct Z2GaugeTheoryDiscretePathIntegralMetropolisParams{F <: AbstractFloat}
    n_side::Int
    n_steps::Int
    β::F
    Δτ::F
    Jxy::F
    Jτ::F
end

function Z2GaugeTheoryDiscretePathIntegralMetropolisParams(J::F, h::F, n_side::Int, n_steps::Int, β::F) where {F <: AbstractFloat}
    Δτ = β / n_steps
    Z2GaugeTheoryDiscretePathIntegralMetropolisParams(n_side, n_steps, β, Δτ, Δτ * J, atanh.(exp.(- 2 * Δτ * h)))
end

plaquatte_op(z2field::DiscretePathIntegralZ2GaugeFieldPeriodicSquare2D, i::Int, τ::Int) = 
    prod(map(x -> z2field[x, τ], plaquatte(z2field.lattice, i)))

flux_average(z2field::DiscretePathIntegralZ2GaugeFieldPeriodicSquare2D, τ::Int) = 
    mean(map(x -> plaquatte_op(z2field, x, τ), sites(z2field.lattice)))

function Δ_plaquatte_term(z2field::DiscretePathIntegralZ2GaugeFieldPeriodicSquare2D, b::Int, τ::Int)
    lattice = z2field.lattice
    - 2 * sum(map(p -> prod(map(x -> z2field[x, τ], p)), plaquatte_shared_bond(lattice, b)))
end

function Δ_temporal_correlation_term(z2field::DiscretePathIntegralZ2GaugeFieldPeriodicSquare2D, b::Int, τ::Int)
    n_steps = z2field.n_steps
    - 2 * z2field[b, τ] * (z2field[b, back_into_range(τ + 1, n_steps)] + z2field[b, back_into_range(τ - 1, n_steps)])
end

"""
Metropolis accept rate.
"""
function accept_rate(z2field::DiscretePathIntegralZ2GaugeFieldPeriodicSquare2D, params::Z2GaugeTheoryDiscretePathIntegralMetropolisParams, b::Int, τ::Int)
    Jxy = params.Jxy
    Jτ = params.Jτ
    exp(Jxy * Δ_plaquatte_term(z2field, b, τ) + Jτ * Δ_temporal_correlation_term(z2field, b, τ))
end

function sweep!(z2field::DiscretePathIntegralZ2GaugeFieldPeriodicSquare2D, params::Z2GaugeTheoryDiscretePathIntegralMetropolisParams, n_sweep::Int; observe = nothing, observe_type::DataType = Any)
    results = observe_type[]
    lattice_bonds = bonds(z2field.lattice)

    for _ in 1 : n_sweep
        for τ in 1 : params.n_steps
            for l in lattice_bonds
                if rand() < accept_rate(z2field, params, l, τ)
                    z2field[l, τ] *= -1 
                end
            end
        end
    end
    if observe !== nothing
        push!(results, observe(z2field))
    end

    results
end

#endregion