#region abstract definitions.

abstract type AbstractLattice{SiteType, CoorType} end

sites(lattice::L) where {S, C, L <: AbstractLattice{S, C}} = @error "Site list not defined for lattice type $L."

site_to_coord(lattice::L, site::S) where {S, C, L <: AbstractLattice{S, C}} = @error "Site to coordinate rules not defined for lattice type $L."

coord_to_site(lattice::L, coord::C) where {S, C, L <: AbstractLattice{S, C}} = @error "Coordinate to site rules not defined for lattice type $L."

function nearest_neighbors(lattice::L, site::S) where {S, C, L <: AbstractLattice{S, C}}
    @error "Nearest neighbors not defined for lattice type $L."
end

function next_nearest_neighbors(lattice::L, site::S) where {S, C, L <: AbstractLattice{S, C}}
    @error "Next nearest neighbors not defined for lattice type $L."
end

function neighbors(lattice::L, site::S, level::I) where {S, C, L <: AbstractLattice{S, C}, I <: Integer}
    if level == 1
        return nearest_neighbors(lattice, site)
    end
    if level == 2
        return next_nearest_neighbors(lattice, site)
    end
    @error "Level $level neighbors not defined for lattice type $L."
end

bond(lattice::L, i::S, j::S) where {S, C, L <: AbstractLattice{S, C}} = @error "Bond between $i and $j is not defined for lattice type $L."

bonds(lattice::L) where {S, C, L <: AbstractLattice{S, C}} = @error "Bond list not defined for lattice type $L."

plaquatte(lattice::L, site::S) where {S, C, L <: AbstractLattice{S, C}} = @error "Plaquatte on $site not defined for lattice type $L."

#endregion

#region Square lattices.

#region vanilla square lattice

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

sites(lattice::L) where {L <: SquareLattice2D} = 1 : lattice.n_side^2

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

#endregion

#region centered square lattice, containing two sublattices.

struct CenteredSquareLattice2D{L <: Integer, SubIdx} <: AbstractLattice{L, Tuple{L, L, SubIdx}}
    n_side::L
    site_list::Vector{Tuple{L, L, SubIdx}}
    inverse_list::Array{L, 3}
end

function CenteredSquareLattice2D(::Type{SubIdx}, n_side::L) where {L <: Integer, SubIdx}
    site_num = n_side * n_side
    site_list = Vector{Tuple{L, L, SubIdx}}(undef, site_num * 2)
    inverse_list = Array{L, 3}(undef, n_side, n_side, 2)
    for x in 1 : n_side
        for y in 1 : n_side 
            site = coord_to_site(n_side, x, y)
            site_list[site] = (x, y, SubIdx(1))
            site_list[site + site_num] = (x, y, SubIdx(2))
            inverse_list[x, y, 1] = site
            inverse_list[x, y, 2] = site + site_num
        end
    end
    CenteredSquareLattice2D{L, SubIdx}(n_side, site_list, inverse_list)
end

@enum CenteredSquareLattice2DSublatticeIndex A = 1 B = 2
function CenteredSquareLattice2D(n_side::L) where {L <: Integer}
    CenteredSquareLattice2D(CenteredSquareLattice2DSublatticeIndex, n_side)
end

function coord_to_site(lattice::CenteredSquareLattice2D{L, SubIdx}, coord::Tuple{L, L, SubIdx})::L where {L <: Integer, SubIdx}
    x, y, sub = coord
    sub = Int(sub)
    lattice.inverse_list[x, y, sub]
end

function site_to_coord(lattice::CenteredSquareLattice2D{L, SubIdx}, site::L)::Tuple{L, L, SubIdx} where {L <: Integer, SubIdx}
    lattice.site_list[site]
end

sites(lattice::CenteredSquareLattice2D) = 1 : 2 * lattice.n_side^2

#endregion

#region periodic square lattice, with nearest neighbors defined.

struct PeriodicSquareLattice2D{L <: Integer} <: AbstractLattice{L, Tuple{L, L}}
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

function nearest_neighbors(n_side::L, i::Tuple{L, L})::Tuple{L, L, L, L} where {L <: Integer}
    x, y = i
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

#endregion

#region square lattice with bonds defined as a dual centered square lattice.

function bond_to_dual_site(n_side::L, i::Tuple{L, L}, j::Tuple{L, L}, ::Type{SubIdx})::Tuple{L, L, SubIdx} where {L <: Integer, SubIdx}
    nn = nearest_neighbors(n_side, i)
    j_site = coord_to_site(n_side, j...)
    if ! (j_site in nn)
        @error "No bond between $i and $j."
    end
    if j_site == nn[1]
        return (i..., SubIdx(2))
    end
    if j_site == nn[2]
        return (j..., SubIdx(2))
    end
    if j_site == nn[3]
        return (i..., SubIdx(1))
    end
    if j_site == nn[4]
        return (j..., SubIdx(1))
    end
end

function bond_to_dual_site(n_side::L, i::Tuple{L, L}, j::Tuple{L, L})::Tuple{L, L} where {L <: Integer}
    bond_to_dual_site(n_side, i, j, CenteredSquareLattice2DSublatticeIndex)    
end

struct PeriodicSquareLattice2DWithBonds{L <: Integer, SubIdx} <: AbstractLattice{L, Tuple{L, L, SubIdx}}
    lattice::PeriodicSquareLattice2D{L}
    bonds_dual_lattice::CenteredSquareLattice2D{L, SubIdx}
    bond_to_dual_site::Matrix{L}
end

function PeriodicSquareLattice2DWithBonds(::Type{SubIdx}, n_side::L) where {L <: Integer, SubIdx}
    lattice = PeriodicSquareLattice2D(n_side)
    bonds_dual_lattice = CenteredSquareLattice2D(SubIdx, n_side)
    bond_mapping = zeros(L, n_side * n_side, n_side * n_side)
    for i in 1 : n_side * n_side
        for j in nearest_neighbors(lattice, i)
            bond_mapping[i, j] = coord_to_site(bonds_dual_lattice, bond_to_dual_site(n_side, site_to_coord(n_side, i), site_to_coord(n_side, j), SubIdx))
        end
    end
    PeriodicSquareLattice2DWithBonds(lattice, bonds_dual_lattice, bond_mapping)
end

PeriodicSquareLattice2DWithBonds(n_side::L) where {L <: Integer} = PeriodicSquareLattice2DWithBonds(CenteredSquareLattice2DSublatticeIndex, n_side)

bond(lattice::PeriodicSquareLattice2DWithBonds{L, SubIdx}, i::L, j::L) where {L <: Integer, SubIdx} = lattice.bond_to_dual_site[i, j]

bonds(lattice::PeriodicSquareLattice2DWithBonds{L, SubIdx}) where {L <: Integer, SubIdx} = sites(lattice.bonds_dual_lattice)

nearest_neighbors(lattice::PeriodicSquareLattice2DWithBonds{L, SubIdx}, i::L) where {L <: Integer, SubIdx} = nearest_neighbors(lattice.lattice, i)

function plaquatte(lattice::PeriodicSquareLattice2DWithBonds{L, SubIdx}, i::L) where {L <: Integer, SubIdx}
    nn_i = nearest_neighbors(lattice, i)
    v1 = i
    v2, v3 = nn_i[1], nn_i[3]
    v4 = nearest_neighbors(lattice, v2)[3]    
    (bond(lattice, v1, v3), bond(lattice, v1, v2), bond(lattice, v3, v4), bond(lattice, v2, v4))
end

function plaquatte_shared_bond(lattice::PeriodicSquareLattice2DWithBonds{L, SubIdx}, b::L) where {L <: Integer, SubIdx}
    n_side = lattice.lattice.lattice.n_side
    if b > n_side
        
    end
end

#endregion

#endregion