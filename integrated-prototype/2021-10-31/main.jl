using ProgressMeter
using Statistics
include("../lib.jl")

##

n_heat = 2000
n_sweep = 1000
n_bin = 10
n_side = 10
show_progress = true 
immediate_store = true

J = 1.0
Δτ = 0.1

T_steps = 80
T_range = LinRange(0.1, 4.1, T_steps)
h_steps = 80
h_range = LinRange(0.2, 3.1, h_steps)

# We are going to run the program by `julia main.jl`
#working_directory = "D:\\Projects\\thermal-z2\\integrated-prototype\\"
working_directory = "./"
output_text_name = "output"

open(working_directory * output_text_name, "w") do file
    println(file, "n_side   =   $n_side")
    println(file, "dtau     =   $Δτ")
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
        n_steps = Int(round(β / Δτ))

        lattice = SquareLattice2DPeriodic(n_side) 
        σ = ones_Z2_gauge_field_DPI(Int, lattice, n_steps)
        model = IsingGaugeTheoryDPIMetropolisMC(Float64, σ, J, h, Δτ)
        
        sweep!(model, σ, n_heat)

        bin = []
        for _ in 1 : n_bin
            push!(bin, mean(sweep!(model, σ, n_sweep, observe = x -> flux_average(x, 1), observable_type = Float64)))
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
