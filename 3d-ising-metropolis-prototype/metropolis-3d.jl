struct AnisotropicIsingSimParams3DMetropolis{N <: Integer}
    nx::N
    ny::N
    nz::N

    n_sweep::N
end

function init_field_3d(
    ::Type{S}, 
    params::AnisotropicIsingSimParams3DMetropolis{N}
)::IsingField3D{S} where {N <: Integer, S <: Integer}
    nx = params.nx
    ny = params.ny
    nz = params.nz

    rand((zero(S), one(S)), (nx, ny, nz))
end

function accept_rate(field::IsingField3D{S}, 
    model::AnisotropicIsingModel3D{F},
    params::AnisotropicIsingSimParams3DMetropolis{N},
    x::N, y::N, z::N
) where {F <: AbstractFloat, S <: Integer, N <: Integer}
    σ = field
    @unpack Jx, Jy, Jz = model
    @unpack nx, ny, nz, n_sweep = params
    
    i = (x, y, z)
    ipx = (back_into_range(x + 1, nx), y, z)
    imx = (back_into_range(x - 1, nx), y, z)
    ipy = (x, back_into_range(y + 1, ny), z)
    imy = (x, back_into_range(y - 1, ny), z)
    ipz = (x, y, back_into_range(z + 1, nz))
    imz = (x, y, back_into_range(z - 1, nz))

    ΔF = (
        2 * Jx * σ[i...] * (σ[ipx...] + σ[imx...]) 
        + 2 * Jy * σ[i...] * (σ[ipy...] + σ[imy...]) 
        + 2 * Jz * σ[i...] * (σ[ipz...] + σ[imz...])
    )

    exp(- ΔF)
end

function run!(field::IsingField3D{S}, 
    model::AnisotropicIsingModel3D{F}, 
    params::AnisotropicIsingSimParams3DMetropolis{N};
    observe = nothing
) where {F <: AbstractFloat, S <: Integer, N <: Integer}
    @unpack nx, ny, nz, n_sweep = params

    sweep_data = []

    for _ in 1 : n_sweep
        for x in 1 : nx
            for y in 1 : ny
                for z in 1 : nz
                    if rand() < accept_rate(field, model, params, x, y, z)
                        field[x, y, z] *= -1
                    end
                end
            end
        end

        if observe !== nothing
            push!(sweep_data, observe(field))
        end
    end

    sweep_data
end