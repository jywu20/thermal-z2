# Note by the author: the architecture in this file now seems awkward and clumsy.
# My mind when I was writing down these codes was soaked in classical OOP patterns, and 
# the relationships between CenteredSquareLattice2D, SquareLattice2D and so on are quite
# like the inheritance hierarchy. 
# This is actually not the best way to do things in Julia.
# A much better way is to use closures to store things like site list, bond list, etc.

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

bond_lattice(lattice::L) where {S, C, L <: AbstractLattice{S, C}} = @error "The lattice constructed by bonds not defined for lattice type $L."

plaquatte(lattice::L, site::S) where {S, C, L <: AbstractLattice{S, C}} = @error "Plaquatte on $site not defined for lattice type $L."

plaquattes(lattice::L) where {S, C, L <: AbstractLattice{S, C}} = @error "Plaquattes not defined for lattice type $L."

n_side(lattice::L) where L <: AbstractLattice = @error "No well defined side length for lattice type $L."

#endregion

#region Square lattices.

#region vanilla square lattice

struct SquareLattice2D <: AbstractLattice{Int, Tuple{Int, Int}}
    n_side::Int
    site_list::Vector{Tuple{Int, Int}}
    inverse_list::Matrix{Int}
end

function SquareLattice2D(n_side::Int)::SquareLattice2D 
    site_list = Vector{Tuple{Int, Int}}(undef, n_side * n_side)
    inverse_list = Matrix{Int}(undef, n_side, n_side)
    for x in 1 : n_side
        for y in 1 : n_side 
            site = coord_to_site(n_side, x, y)
            site_list[site] = (x, y)
            inverse_list[x, y] = site
        end
    end
    SquareLattice2D(n_side, site_list, inverse_list)
end

n_side(lattice::SquareLattice2D) = lattice.n_side

sites(lattice::SquareLattice2D) = 1 : lattice.n_side^2

function site_to_coord(n_side::Int, i::Int)::Tuple{Int, Int}
    y = i % n_side
    if y == 0
        y = n_side
    end
    x = 1 + Int((i - y) / n_side)
    
    (x, y)
end

function coord_to_site(n_side::Int, x::Int, y::Int)::Int
    (x - 1) * n_side + y
end

function site_to_coord(lattice::SquareLattice2D, site::Int)::Tuple{Int, Int}
    lattice.site_list[site]
end

function coord_to_site(lattice::SquareLattice2D, site::Tuple{Int, Int})::Int
    lattice.inverse_list[site...]
end

#endregion

#region centered square lattice, containing two sublattices.

@enum CenteredSquareLattice2DSublatticeIndex A = 1 B = 2

struct CenteredSquareLattice2D <: AbstractLattice{Int, Tuple{Int, Int, CenteredSquareLattice2DSublatticeIndex}}
    n_side::Int
    site_list::Vector{Tuple{Int, Int, CenteredSquareLattice2DSublatticeIndex}}
    inverse_list::Array{Int, 3}
end

function CenteredSquareLattice2D(n_side::Int)
    site_num = n_side * n_side
    site_list = Vector{Tuple{Int, Int, CenteredSquareLattice2DSublatticeIndex}}(undef, site_num * 2)
    inverse_list = Array{Int, 3}(undef, n_side, n_side, 2)
    for x in 1 : n_side
        for y in 1 : n_side 
            site = coord_to_site(n_side, x, y)
            site_list[site] = (x, y, A)
            site_list[site + site_num] = (x, y, B)
            inverse_list[x, y, 1] = site
            inverse_list[x, y, 2] = site + site_num
        end
    end
    CenteredSquareLattice2D(n_side, site_list, inverse_list)
end

function coord_to_site(lattice::CenteredSquareLattice2D, coord::Tuple{Int, Int, CenteredSquareLattice2DSublatticeIndex})::Int
    x, y, sub = coord
    sub = Int(sub)
    lattice.inverse_list[x, y, sub]
end

function site_to_coord(lattice::CenteredSquareLattice2D, site::Int)::Tuple{Int, Int, CenteredSquareLattice2DSublatticeIndex}
    lattice.site_list[site]
end

sites(lattice::CenteredSquareLattice2D) = 1 : 2 * lattice.n_side^2

#endregion

#region periodic square lattice, with nearest neighbors defined.

struct PeriodicSquareLattice2D <: AbstractLattice{Int, Tuple{Int, Int}}
    lattice::SquareLattice2D
    nn_list::Vector{Tuple{Int, Int, Int, Int}}
end

function PeriodicSquareLattice2D(n_side::Int)::PeriodicSquareLattice2D
    lattice = SquareLattice2D(n_side)
    nn_list = Vector{Tuple{Int, Int, Int, Int}}(undef, n_side * n_side)
    for site in 1 : n_side * n_side
        nn_list[site] = nearest_neighbors(n_side, site)
    end
    PeriodicSquareLattice2D(lattice, nn_list)
end

function nearest_neighbors(n_side::Int, i::Int)::Tuple{Int, Int, Int, Int}
    x, y = site_to_coord(n_side, i)
    zero_to_max(i) = i > 0 ? i : i + n_side
    (
        coord_to_site(n_side, zero_to_max(mod(x + 1, n_side)), y),
        coord_to_site(n_side, zero_to_max(mod(x - 1, n_side)), y),
        coord_to_site(n_side, x, zero_to_max(mod(y + 1, n_side))),
        coord_to_site(n_side, x, zero_to_max(mod(y - 1, n_side)))
    )
end

function nearest_neighbors(n_side::Int, i::Tuple{Int, Int})::Tuple{Int, Int, Int, Int}
    x, y = i
    zero_to_max(i) = i > 0 ? i : i + n_side
    (
        coord_to_site(n_side, zero_to_max(mod(x + 1, n_side)), y),
        coord_to_site(n_side, zero_to_max(mod(x - 1, n_side)), y),
        coord_to_site(n_side, x, zero_to_max(mod(y + 1, n_side))),
        coord_to_site(n_side, x, zero_to_max(mod(y - 1, n_side)))
    )
end

function nearest_neighbors(lattice::PeriodicSquareLattice2D, site::Int)::Tuple{Int, Int, Int, Int} 
    lattice.nn_list[site]
end

site_x_forward_move(lattice::PeriodicSquareLattice2D, site::Int)::Int = nearest_neighbors(lattice, site)[1]
site_x_backward_move(lattice::PeriodicSquareLattice2D, site::Int)::Int = nearest_neighbors(lattice, site)[2]
site_y_forward_move(lattice::PeriodicSquareLattice2D, site::Int)::Int = nearest_neighbors(lattice, site)[3]
site_y_backward_move(lattice::PeriodicSquareLattice2D, site::Int)::Int = nearest_neighbors(lattice, site)[4]

sites(lattice::PeriodicSquareLattice2D) = sites(lattice.lattice)

n_side(lattice::PeriodicSquareLattice2D) = n_side(lattice.lattice)

#endregion

#region square lattice with bonds defined as a dual centered square lattice.

function bond_to_dual_site(n_side::Int, i::Tuple{Int, Int}, j::Tuple{Int, Int})::Tuple{Int, Int, CenteredSquareLattice2DSublatticeIndex} 
    nn = nearest_neighbors(n_side, i)
    j_site = coord_to_site(n_side, j...)
    if ! (j_site in nn)
        @error "No bond between $i and $j."
    end
    if j_site == nn[1]
        return (i..., B)
    end
    if j_site == nn[2]
        return (j..., B)
    end
    if j_site == nn[3]
        return (i..., A)
    end
    if j_site == nn[4]
        return (j..., A)
    end
end

struct PeriodicSquareLattice2DWithBonds <: AbstractLattice{Int, Tuple{Int, Int}}
    lattice::PeriodicSquareLattice2D
    bonds_dual_lattice::CenteredSquareLattice2D
    bond_to_dual_site::Matrix{Int}
end

function PeriodicSquareLattice2DWithBonds(n_side::Int)
    lattice = PeriodicSquareLattice2D(n_side)
    bonds_dual_lattice = CenteredSquareLattice2D(n_side)
    bond_mapping = zeros(Int, n_side * n_side, n_side * n_side)
    for i in 1 : n_side * n_side
        for j in nearest_neighbors(lattice, i)
            bond_mapping[i, j] = coord_to_site(bonds_dual_lattice, bond_to_dual_site(n_side, site_to_coord(n_side, i), site_to_coord(n_side, j)))
        end
    end
    PeriodicSquareLattice2DWithBonds(lattice, bonds_dual_lattice, bond_mapping)
end

sites(lattice::PeriodicSquareLattice2DWithBonds) = sites(lattice.lattice)

n_side(lattice::PeriodicSquareLattice2DWithBonds) = n_side(lattice.lattice)

bond(lattice::PeriodicSquareLattice2DWithBonds, i::Int, j::Int) = lattice.bond_to_dual_site[i, j]

bonds(lattice::PeriodicSquareLattice2DWithBonds) = sites(lattice.bonds_dual_lattice)

bond_lattice(lattice::PeriodicSquareLattice2DWithBonds) = lattice.bonds_dual_lattice

nearest_neighbors(lattice::PeriodicSquareLattice2DWithBonds, i::Int) = nearest_neighbors(lattice.lattice, i)

site_x_forward_move(lattice::PeriodicSquareLattice2DWithBonds, site::Int)::Int = site_x_forward_move(lattice.lattice, site)
site_x_backward_move(lattice::PeriodicSquareLattice2DWithBonds, site::Int)::Int = site_x_backward_move(lattice.lattice, site)
site_y_forward_move(lattice::PeriodicSquareLattice2DWithBonds, site::Int)::Int = site_y_forward_move(lattice.lattice, site)
site_y_backward_move(lattice::PeriodicSquareLattice2DWithBonds, site::Int)::Int = site_y_backward_move(lattice.lattice, site)

function plaquattes(lattice::PeriodicSquareLattice2DWithBonds)
    sites(lattice)
end

function plaquatte(lattice::PeriodicSquareLattice2DWithBonds, i::Int)
    nn_i = nearest_neighbors(lattice, i)
    v1 = i
    v2, v3 = nn_i[1], nn_i[3]
    v4 = nearest_neighbors(lattice, v2)[3]    
    (bond(lattice, v1, v3), bond(lattice, v1, v2), bond(lattice, v3, v4), bond(lattice, v2, v4))
end

function plaquatte_shared_bond(lattice::PeriodicSquareLattice2DWithBonds, b::Int)
    n_sites = length(sites(lattice))
    # B-bonds, or bonds in the x direction
    if b > n_sites
        top_point = b - n_sites
        left_top_point = site_y_backward_move(lattice, top_point)
        return (plaquatte(lattice, left_top_point), plaquatte(lattice, top_point))
    end
    # A-bonds or bonds in the y direction
    left_point = b
    left_top_point = site_x_backward_move(lattice, b)
    return (plaquatte(lattice, left_top_point), plaquatte(lattice, left_point))
end

#endregion

#endregion