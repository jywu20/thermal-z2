struct DiscretePathIntegralZ2GaugeFieldConfiguration{L <: AbstractLattice, D, N}
    lattice::L
    data::Array{D, N}
end

getindex(z2field::DiscretePathIntegralZ2GaugeFieldConfiguration, idx...) = z2field.data[idx...]
setindex!(z2field::DiscretePathIntegralZ2GaugeFieldConfiguration, value, idx) = z2field.data[idx...] = value

const DiscretePathIntegralZ2GaugeFieldConfiguration2D{L <: AbstractLattice, D} = DiscretePathIntegralZ2GaugeFieldConfiguration{L, D, 3}

function DiscretePathIntegralZ2GaugeFieldConfiguration2D(::Type{D}, lattice::L, n_steps::L) where {L <: AbstractLattice, D} 
    n_side = length(sites(lattice))
    DiscretePathIntegralZ2GaugeFieldConfiguration2D(lattice, Array{D, 3}(undef, n_side, n_side, n_steps))
end

function flip!(z2field::DiscretePathIntegralZ2GaugeFieldConfiguration2D{L, D}, i::S) where {S, C, D, L <: AbstractLattice{S, C}}
    z2field[i] *= -1
end

struct Z2GaugeFieldModel{F <: AbstractFloat}
    J::F
    h::F
end

struct Z2GaugeFieldDiscretePathIntegralParams{N <: Integer, F <: AbstractFloat}
    n_side::N
    n_steps::N
    β::F
end

function accept_rate(z2field::DiscretePathIntegralZ2GaugeFieldConfiguration2D, i)
    z2field[i]
end