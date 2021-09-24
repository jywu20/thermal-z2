abstract type AbstractLattice{SiteType, CoorType} end

sites(lattice::L) where {S, C, L <: AbstractLattice{S, C}} = @error "Site list not defined for lattice type $L."

site_to_coord(lattice::L, site::S) where {S, C, L <: AbstractLattice{S, C}} = @error "Site to coordinate rules not defined for lattice type $L."

coord_to_site(lattice::L, coord::C) where {S, C, L <: AbstractLattice{S, C}} = @error "Coordinate to site rules not defined for lattice type $L."

struct SquareLattice2D{L <: Integer} <: AbstractLattice{L, Tuple{L, L}}
    n_side::L
    site_list::Vector{Tuple{L, L}}
    inverse_list::Matrix{L}
end

function SquareLattice2D(n_side::L)::SquareLattice2D{L} where {L <: Integer}
    site_list = Vector{Tuple{L, L}}(undef, n_side * n_side)
    inverse_list = Matrix{L}(undef, n_side, n_side)
    for x in 1 : n_side
        for y in 1 : n_side 
            site = coord_to_site(n_side, x, y)
            site_list[site] = (x, y)
            inverse_list[x, y] = site
        end
    end
    SquareLattice2D{L}(n_side, site_list, inverse_list)
end

function site_to_coord(n_side::L, i::L)::Tuple{L, L} where {L <: Integer}
    y = i % n_side
    if y == 0
        y = n_side
    end
    x = 1 + Int((i - y) / n_side)
    
    (x, y)
end

function coord_to_site(n_side::L, x::L, y::L)::L where {L <: Integer}
    (x - 1) * n_side + y
end

function site_to_coord(lattice::SquareLattice2D{L}, site::L)::Tuple{L, L} where {L <: Integer}
    lattice.site_list[site]
end
function coord_to_site(lattice::SquareLattice2D{L}, site::Tuple{L, L})::L where {L <: Integer}
    lattice.inverse_list[site...]
end

abstract type AbstractLatticeWithNeighbors{SiteType, CoorType, NNType} <: AbstractLattice{SiteType, CoorType} end

function nearest_neighbors(lattice::L, site::S)::N where {S, C, N, L <: AbstractLatticeWithNeighbors{S, C, N}}
    @error "Nearest neighbors not defined for lattice type $L"
end

function next_nearest_neighbors(lattice::L, site::S)::N where {S, C, N, L <: AbstractLatticeWithNeighbors{S, C, N}}
    @error "Next nearest neighbors not defined for lattice type $L"
end

function neighbors(lattice::L, site::S, level::I)::N where {S, C, N, L <: AbstractLatticeWithNeighbors{S, C, N}, I <: Integer}
    if level == 1
        return nearest_neighbors(lattice, site)
    end
    if level == 2
        return next_nearest_neighbors(lattice, site)
    end
    @error "Level $level neighbors not defined for lattice type $L"
end

struct PeriodicSquareLattice2D{L <: Integer} <: AbstractLatticeWithNeighbors{L, Tuple{L, L}, Tuple{L, L, L, L}}
    lattice::SquareLattice2D{L}
    nn_list::Vector{Tuple{L, L, L, L}}
end

function PeriodicSquareLattice2D(n_side::L)::PeriodicSquareLattice2D{L} where {L <: Integer}
    lattice = SquareLattice2D(n_side)
    nn_list = Vector{Tuple{L, L, L, L}}(undef, n_side * n_side)
    for site in 1 : n_side * n_side
        nn_list[site] = nearest_neighbors(n_side, site)
    end
    PeriodicSquareLattice2D(lattice, nn_list)
end

function nearest_neighbors(n_side::L, i::L)::Tuple{L, L, L, L} where {L <: Integer}
    x, y = site_to_coord(n_side, i)
    zero_to_max(i) = i > 0 ? i : i + n_side
    (
        coord_to_site(n_side, zero_to_max(mod(x + 1, n_side)), y),
        coord_to_site(n_side, zero_to_max(mod(x - 1, n_side)), y),
        coord_to_site(n_side, x, zero_to_max(mod(y + 1, n_side))),
        coord_to_site(n_side, x, zero_to_max(mod(y - 1, n_side)))
    )
end

function nearest_neighbors(lattice::PeriodicSquareLattice2D{L}, site::L)::Tuple{L, L, L, L} where {L <: Integer}
    lattice.nn_list[site]
end
