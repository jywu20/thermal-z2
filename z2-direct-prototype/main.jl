using ProgressMeter
using Statistics
include("z2.jl")

##

n_heat = 1000
n_sweep = 1000
n_bin = 10
n_steps = 20
n_side = 10
show_progress = true
immediate_store = true

J = 1.0

T_steps = 60
T_range = LinRange(0.1, 10.1, T_steps)
h_steps = 50
h_range = LinRange(0.1, 4.1, h_steps)
Δh = (h_range[end] - h_range[1]) / h_steps

working_directory = "D:\\Projects\\thermal-z2\\z2-direct-prototype\\"
output_text_name = "output-2021-10-8"

open(working_directory * output_text_name, "w") do file
    println(file, "n_side   =   $n_side")
    println(file, "n_steps  =   $n_steps")
    println(file, "n_bin    =   $n_bin")
    println(file, "n_heat   =   $n_heat")
    println(file, "n_sweep  =   $n_sweep")
    println(file, "T_min    = $(T_range[1]),    T_max = $(T_range[end]),   steps = $T_steps")
    println(file, "h_min    = $(h_range[1]),    T_max = $(h_range[end]),   steps = $h_steps")
    println(file, "J        =   $J")
end

flux_history = zeros(T_steps, h_steps)

progress = Progress(T_steps * h_steps)

for (i, T) in enumerate(T_range)
    for (j, h) in enumerate(h_range)
        β = 1 / T
        lattice = PeriodicSquareLattice2DWithBonds(n_side)
        z2field = all_ones_discrete_path_integral_z2_gauge_field(lattice, n_steps)
        params = Z2GaugeTheoryDiscretePathIntegralMetropolisParams(J, h, n_side, n_steps, β)
        
        sweep!(z2field, params, n_heat)

        bin = []
        for _ in 1 : n_bin
            push!(bin, mean(sweep!(z2field, params, n_sweep, observe = x -> flux_average(x, 1), observe_type = Float64)))
        end
        flux_history[i, j] = mean(bin)

        if immediate_store
            open(working_directory * output_text_name, "a") do file
                println(file, "T = $T       h = $h      mean_flux = $(flux_history[i, j])")
            end
        end

        if show_progress
            next!(progress)
        end
    end
end
