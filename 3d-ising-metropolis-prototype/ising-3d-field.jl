# In this project only square lattice is investigated, therefore we do not put the word "lattice"
# into type names.

"""
A 3d classical Ising field.

We use the convention that the first index represents the x coordinate, etc.
"""
const IsingField3D = Array{T, 3} where {T <: Integer}

"""
A 2d transverse field Ising field is just a 3d classical Ising field, because the field configuration 
of the path integral of the former one is just the latter.

We use the convention that the first index represents the time, the second index the x coordinate and
the third index the y coordinate.
"""
const TFIsingField2D = IsingField3D

struct AnisotropicIsingModel3D{F <: AbstractFloat}
    Jx::F
    Jy::F
    Jz::F
end

struct TransverseFieldIsingModel2D{F <: AbstractFloat}
    J::F
    h::F
end

"""
Map the parameters of a 2d transverse field Ising model (TFIM) into ones of a 3d classical Ising model (CIM).
"""
function tfim_2d_to_cim_3d(model::TransverseFieldIsingModel2D{F}) where {F <: AbstractFloat}
    
end