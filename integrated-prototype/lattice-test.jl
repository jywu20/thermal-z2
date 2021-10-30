include("lib.jl")

##

n_side = 4
square_lattice = SquareLattice2DPeriodic(n_side)

##
 
site_to_coord(square_lattice, 13)

##

coord_to_site(square_lattice, (4, 3))

##
bond_to_sites(square_lattice, 28)

##

sites_to_bond(square_lattice, 3, 4)

##
nearest_neighbors(square_lattice, 1)

##
bonds_around_site(square_lattice, 8)

##

plaquatte_to_corners(square_lattice, 12)

##
plaquatte_to_bonds(square_lattice, 16)

##
plaquatte_containing_bond(square_lattice, 27)