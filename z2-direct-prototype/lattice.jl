abstract type AbstractLattice{SiteType, CoorType} end

struct SquareLattice2D{L <: Integer} <: AbstractLattice{L, Tuple{L, L}}
    n_side::L
    site_list::Vector{Tuple{L, L}}
    inverse_list::Matrix{L}
end

