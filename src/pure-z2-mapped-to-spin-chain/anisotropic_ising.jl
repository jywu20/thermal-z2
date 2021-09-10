struct AnisotropicIsing2DHamParams 
    kx::Float64
    ky::Float64
end

struct AnisotropicIsing2DSimParams 
    nx::Int64
    ny::Int64
    n_sweep::Int64
end

const IsingField2D = Matrix{Int64}

function init(sim_params::AnisotropicIsing2DSimParams)
    nx = sim_params.nx
    ny = sim_params.ny
    ising_field = rand((-1, 1), (nx, ny))
    ising_field
end

function flip!(ising_field::IsingField2D, cluster)
    for pos in cluster
        ising_field[pos...] *= -1
    end
end

function run!(ham_params::AnisotropicIsing2DHamParams, 
    sim_params::AnisotropicIsing2DSimParams, 
    ising_field::IsingField2D; 
    observe = nothing)::Vector
    
    kx = ham_params.kx
    ky = ham_params.ky

    px = 1.0 - exp(- 2kx)
    py = 1.0 - exp(- 2ky)

    nx = sim_params.nx
    ny = sim_params.ny
    n_sweep = sim_params.n_sweep

    sweep_data = []

    for _ in 1 : n_sweep

        start_point_x = rand(1:nx)
        start_point_y = rand(1:ny)

        pocket = [(start_point_x, start_point_y)]
        cluster = [(start_point_x, start_point_y)]

        while length(pocket) != 0
            this_point = pop!(pocket)
            this_x, this_y = this_point

            next_y = this_y
            for next_x in [back_into_range(this_x + 1, nx), back_into_range(this_x - 1, nx)]
                if ising_field[this_x, this_y] == ising_field[next_x, next_y] && !((next_x, next_y) in cluster)
                    if rand(Float64) < px
                        push!(pocket, (next_x, next_y))
                        push!(cluster, (next_x, next_y))
                    end
                end
            end

            next_x = this_x
            for next_y in [back_into_range(this_y + 1, ny), back_into_range(this_y - 1, ny)]
                if ising_field[this_x, this_y] == ising_field[next_x, next_y] && !((next_x, next_y) in cluster)
                    if rand(Float64) < py
                        push!(pocket, (next_x, next_y))
                        push!(cluster, (next_x, next_y))
                    end
                end
            end

        end

        flip!(ising_field, cluster)

        if observe !== nothing
            push!(sweep_data, observe(ising_field))
        end
    end

    sweep_data
end