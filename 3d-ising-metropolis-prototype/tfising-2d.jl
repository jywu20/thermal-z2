"""
A 2d transverse field Ising field is just a 3d classical Ising field, because the field configuration 
of the path integral of the former one is just the latter.

We use the convention that the last index represents the time, the first index the x coordinate and
the second index the y coordinate. Or in other words, the z coordinate is actually τ.
"""
const TFIsingFieldWolff2D = IsingField3D

struct TransverseFieldIsingModel2D{F <: AbstractFloat}
    J::F
    h::F
end

struct TransverseFieldIsingSimParams2D{F <: AbstractFloat, N <: Integer}
    β::F
    n_steps::N
end

"""
Map the parameters of a 2d transverse field Ising model (TFIM) into ones of a 3d classical Ising model (CIM).
"""
function tfim_2d_to_cim_3d(
    model::TransverseFieldIsingModel2D{F}, 
    sim_params::TransverseFieldIsingSimParams2D{F, N}
)::AnisotropicIsingModel3D{F} where {F <: AbstractFloat, N <: Integer}
    β = sim_params.β
    n_steps = sim_params.n_steps
    Δτ = β / n_steps
    J = model.J
    h = model.h

    Jxy = Δτ * J
    Jτ = - log(tanh(Δτ * h)) / 2

    AnisotropicIsingModel3D(Jxy, Jxy, Jτ)
end