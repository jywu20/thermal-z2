function tight_binding_kinetic_matrix_periodic(lattice::PeriodicSquareLattice2DWithBonds)
    lattice_sites = sites(lattice)
    T = zeros(length(lattice_sites), length(lattice_sites))
    for i in lattice_sites
        for j in nearest_neighbors(lattice, i)
            T[i, j] = 1.0
        end
    end
    T
end

function z2_to_tight_binding_kinetic_matrix(z2field::DiscretePathIntegralZ2GaugeFieldPeriodicSquare2D, τ::Int)
    lattice_sites = sites(z2field.lattice)
    T = zeros(length(lattice_sites), length(lattice_sites))
    for i in lattice_sites
        for j in nearest_neighbors(lattice, i) 
            T[i, j] = z2field[i, j, τ]
        end
    end
    T
end

function Tij(n_sites, i, j)
    res = zeros(n_sites, n_sites)
    res[i, j] = 1
    res[j, i] = 1
    res
end

"""
Returns plaquattes that are included in the checkboard decomposition.
The lattice should have an even n_side.
"""
function checkboard_decomposition(lattice::PeriodicSquareLattice2DWithBonds)
    n_side_lattice = n_side(lattice)
    sub_A = Int[i + j for i in 0 : 2n_side_lattice : n_side_lattice^2 - 1 for j in 1 : 2 : n_side_lattice]
    sub_B = Int[i + j for i in n_side_lattice : 2n_side_lattice : n_side_lattice^2 for j in 2 : 2 : n_side_lattice]
    (sub_A, sub_B)
end

function plaquatte_hamiltonian(lattice::PeriodicSquareLattice2DWithBonds, p::Int)
    plaquatte(lattice, p)    # TODO
end

function checkboard_decomposition_hamiltonians(lattice::PeriodicSquareLattice2DWithBonds)
    # TODO
end