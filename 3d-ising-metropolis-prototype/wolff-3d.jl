# In this project only square lattice is investigated, therefore we do not put the word "lattice"
# into type names.

"""
A 3d classical Ising field.

We use the convention that the first index represents the x coordinate, etc.
"""
const IsingField3D = Array{T, 3} where {T <: Integer}

struct AnisotropicIsingModel3D{F <: AbstractFloat}
    Jx::F
    Jy::F
    Jz::F
end

struct AnisotropicIsingSimParams3DWolff{N <: Integer}
    nx::N
    ny::N
    nz::N

    n_sweep::N
end

function init_field_3d(
    ::Type{S}, 
    params::AnisotropicIsingSimParams3DWolff{N}
)::IsingField3D{S} where {N <: Integer, S <: Integer}
    nx = params.nx
    ny = params.ny
    nz = params.nz

    rand((zero(S), one(S)), (nx, ny, nz))
end

function flip!(ising_field::IsingField3D, cluster)
    for pos in cluster
        ising_field[pos...] *= -1
    end
end

function run!(
    field::IsingField3D{S}, 
    model::AnisotropicIsingModel3D{F}, 
    params::AnisotropicIsingSimParams3DWolff{N};
    observe = nothing
) where {F <: AbstractFloat, S <: Integer, N <: Integer}
    Jx = model.Jx
    Jy = model.Jy
    Jz = model.Jz

    nx = params.nx
    ny = params.ny
    nz = params.nz

    n_sweep = params.n_sweep

    px = 1.0 - exp(- 2Jx)
    py = 1.0 - exp(- 2Jy)
    pz = 1.0 - exp(- 2Jz)

    sweep_data = []

    for _ in 1 : n_sweep
        start_point_x = rand(1:nx)
        start_point_y = rand(1:ny)
        start_point_z = rand(1:nz)

        pocket = [(start_point_x, start_point_y, start_point_z)]
        cluster = [(start_point_x, start_point_y, start_point_z)]

        while length(pocket) != 0
            this_point = pop!(pocket)
            this_x, this_y, this_z = this_point

            next_y = this_y
            next_z = this_z
            for next_x in [back_into_range(this_x + 1, nx), back_into_range(this_x - 1, nx)]
                if field[this_x, this_y, this_z] == field[next_x, next_y, next_z] && !((next_x, next_y, next_z) in cluster)
                    if rand(Float64) < px
                        push!(pocket, (next_x, next_y, next_z))
                        push!(cluster, (next_x, next_y, next_z))
                    end
                end
            end

            next_x = this_x
            next_z = this_z
            for next_y in [back_into_range(this_y + 1, ny), back_into_range(this_y - 1, ny)]
                if field[this_x, this_y, this_z] == field[next_x, next_y, next_z] && !((next_x, next_y, next_z) in cluster)
                    if rand(Float64) < py
                        push!(pocket, (next_x, next_y, next_z))
                        push!(cluster, (next_x, next_y, next_z))
                    end
                end
            end

            next_x = this_x
            next_y = this_y
            for next_z in [back_into_range(this_z + 1, nz), back_into_range(this_z - 1, nz)]
                if field[this_x, this_y, this_z] == field[next_x, next_y, next_z] && !((next_x, next_y, next_z) in cluster)
                    if rand(Float64) < pz
                        push!(pocket, (next_x, next_y, next_z))
                        push!(cluster, (next_x, next_y, next_z))
                    end
                end
            end
        end

        flip!(field, cluster)

        if observe !== nothing
            push!(sweep_data, observe(field))
        end
    end

    sweep_data
end