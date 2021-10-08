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

struct TFIMDiscretePathIntegralMetropolisParams{F <: AbstractFloat}
    n_side::Int
    n_steps::Int
    β::F
    Δτ::F
    Jxy::F
    Jτ::F
end

function TFIMDiscretePathIntegralMetropolisParams(J::F, h::F, n_side::Int, n_steps::Int, β::F) where {F <: AbstractFloat}
    Δτ = β / n_steps
    TFIMDiscretePathIntegralMetropolisParams(n_side, n_steps, β, Δτ, Δτ * J, atanh.(exp.(- 2 * Δτ * h)))
end

function Δ_spacial_correlation_term(ising_field::DiscretePathIntegralIsingFieldPerdiodicSquare2D, i::Int, τ::Int)
    lattice = ising_field.lattice
    -2 * ising_field[i, τ] * sum(map(x -> ising_field[x, τ], nearest_neighbors(lattice, i)))
end

function Δ_temporal_correlation_term(ising_field::DiscretePathIntegralIsingFieldPerdiodicSquare2D, i::Int, τ::Int)
    n_steps = ising_field.n_steps
    - 2 * ising_field[i, τ] * (ising_field[i, back_into_range(τ + 1, n_steps)] + ising_field[i, back_into_range(τ - 1, n_steps)])
end

"""
Metropolis accept rate.
"""
function accept_rate(
    ising_field::DiscretePathIntegralIsingFieldPerdiodicSquare2D, 
    params::TFIMDiscretePathIntegralMetropolisParams, 
    i::Int, τ::Int
)
    Jxy = params.Jxy
    Jτ = params.Jτ
    exp(Jxy * Δ_spacial_correlation_term(ising_field, i, τ) + Jτ * Δ_temporal_correlation_term(ising_field, i, τ))
end

function sweep!(
    ising_field::DiscretePathIntegralIsingFieldPerdiodicSquare2D,
    params::TFIMDiscretePathIntegralMetropolisParams,
    n_sweep::Int;
    observe = nothing, observe_type::DataType = Any
)
    results = observe_type[]
    lattice_sites = sites(ising_field.lattice)    

    for _ in 1 : n_sweep
        for τ in 1 : params.n_steps
            for l in lattice_sites
                if rand() < accept_rate(ising_field, params, l, τ)
                    ising_field[l, τ] *= -1 
                end
            end
        end
    end
    if observe !== nothing
        push!(results, observe(ising_field))
    end

    results
end

#endregion