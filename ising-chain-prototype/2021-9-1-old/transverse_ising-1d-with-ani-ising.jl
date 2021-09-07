using Random
using Statistics
using ProgressMeter
using Plots
using LaTeXStrings

function back_into_range(idx, upper)
    if idx > upper
        return idx % upper
    end
    (idx - upper) % upper + upper
end

function anisotropic_ising(;kx, ky, nx, ny, n_bin, n_heat, n_sweep, observe)

    px = 1.0 - exp(- 2kx)
    py = 1.0 - exp(- 2ky)

    function flip(cluster)
        for pos in cluster
            ising_field[pos...] *= -1
        end
    end

    function sweep(time, is_observed)
        sweep_data = []

        for sweep_count in 1:time

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

            flip(cluster)

            if is_observed
                push!(sweep_data, observe(ising_field))
            end
        end

        sweep_data
    end

    ising_field = rand((-1, 1), (nx, ny))

    sweep(n_heat, false)

    bin_data = []
    for bin_count in 1:n_bin
        sweep_data = sweep(n_sweep, true)
        push!(bin_data, mean(sweep_data))
    end

    bin_data
end

##

β = 1000
J = 1.0
h = 0.0

nx = 100
ny = 10
n_heat = 10
n_bin = 10
n_sweep = 400
dh = 0.05

T_range = LinRange(0.05, 4, 20)
# T_range = [0.1]
h_range = LinRange(0.01, 2.5, 40)
mag_to_T_h = Array{Float64}(undef, length(T_range), length(h_range))
mag_to_T_h_dh = Array{Float64}(undef, length(T_range), length(h_range))
progress = Progress(length(T_range) * length(h_range))

for (i, T) in enumerate(T_range), (j, h) in enumerate(h_range)
    β = 1.0 / T
    dtau = β / nx
    j_tau = -log(tanh(dtau * h)) / 2
    result = anisotropic_ising(kx=j_tau, ky=J * dtau, nx=nx, ny=ny, n_bin=n_bin, n_heat=n_heat, n_sweep=n_sweep, observe=abs ∘ mean)
    mag_to_T_h[i, j] = mean(result)

    j_tau_dh = -log(tanh(dtau * (h + dh))) / 2
    result = anisotropic_ising(kx=j_tau_dh, ky=J * dtau, nx=nx, ny=ny, n_bin=n_bin, n_heat=n_heat, n_sweep=n_sweep, observe=abs ∘ mean)
    mag_to_T_h_dh[i, j] = mean(result)

    next!(progress)
end

##

heatmap(h_range, T_range, mag_to_T_h, 
    xlabel = L"h", 
    ylabel = L"T",
    title = "Magnetization"
)

##

heatmap(h_range, T_range, (mag_to_T_h_dh - mag_to_T_h) / dh,
    xlabel = L"h", 
    ylabel = L"T",
    title = "Magnetic susceptibility"
)
