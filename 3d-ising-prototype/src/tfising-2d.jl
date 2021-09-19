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

struct TransverseFieldIsingSimParams2DWolff{F <: AbstractFloat, N <: Integer}
    β::F
    n_steps::N
    n_side::N
    n_sweep::N
end

struct TransverseFieldIsingSimParams2DMetropolis{F <: AbstractFloat, N <: Integer}
    β::F
    n_steps::N
    n_side::N
    n_sweep::N
end

const TransverseFieldIsingSimParams2DNaivePathIntegral{F <: AbstractFloat, N <: Integer} = Union{
    TransverseFieldIsingSimParams2DMetropolis{F, N}, 
    TransverseFieldIsingSimParams2DWolff{F, N} 
}

"""
Map the parameters of a 2d transverse field Ising model (TFIM) into ones of a 3d classical Ising model (CIM).
"""
function tfim_2d_to_cim_3d(
    model::TransverseFieldIsingModel2D{F}, 
    sim_params::Sim
)::AnisotropicIsingModel3D{F} where {
    F <: AbstractFloat, N <: Integer, 
    Sim <: TransverseFieldIsingSimParams2DNaivePathIntegral{F, N}
}
    β = sim_params.β
    n_steps = sim_params.n_steps
    Δτ = β / n_steps
    J = model.J
    h = model.h

    Jxy = Δτ * J
    Jτ = - log(tanh(Δτ * h)) / 2

    AnisotropicIsingModel3D(Jxy, Jxy, Jτ)
end

run!(
    field::TFIsingFieldWolff2D{S},
    model::TransverseFieldIsingModel2D{F},
    params::TransverseFieldIsingSimParams2DWolff{F, N};
    observe = nothing
) where {F <: AbstractFloat, S <: Integer, N <: Integer} = 
    run!(field, tfim_2d_to_cim_3d(model, params), 
        AnisotropicIsingSimParams3DWolff(params.n_side, params.n_side, params.n_steps, params.n_sweep), 
        observe = observe
    )