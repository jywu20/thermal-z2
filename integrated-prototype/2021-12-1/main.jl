using Statistics
using ProgressMeter

include("../lib.jl")

###################################################################################
# Parameters
###################################################################################

n_sweep = 100
n_heat = 100
n_bin = 10
n_side = 4
n_τ = 100
n_wrap = 10

Δτ = 0.01
β = n_τ * Δτ 

t = 1.0
h_Z2 = 1.0
h_TFIM = 1.0
J = 1.0
g = 1.0

show_progress = true
working_path = "D:/Projects/thermal-z2/integrated-prototype/2021-12-1/"

###################################################################################
# Lattice and fields initialization.
###################################################################################

lattice = SquareLattice2DPeriodic(n_side)
lattice_sites = sites(lattice)
lattice_bonds = bonds(lattice)

σ = ones_Z2_gauge_field_DPI(Int, lattice, n_τ)
model_Z2 = IsingGaugeTheoryDPIMetropolisMC(Float64, σ, g, h_Z2, Δτ)
model_Z2_fermion = Z2SpinlessFermionSimpleDQMC(Float64, σ, n_wrap, t, Δτ)
aux = Z2SpinlessFermionSimpleAuxField(Float64, model_Z2_fermion, σ)

s = ones_Ising_field_DPI(Int, lattice, n_τ)
model_TFIM = TransverseFieldIsingModelDPIMetropolisMC(Float64, s, 0.0, h_TFIM, Δτ)

model_coupling = Z2IsingMinimalCouplingDPI{Float64}(J, Δτ)

###################################################################################
# Parameters refreshing
###################################################################################

function refresh_params(t′, h_Z2′, h_TFIM′, J′, g′, T′)
    global t = t′
    global h_Z2 = h_Z2′
    global h_TFIM = h_TFIM′
    global J = J′
    global g = g′
    n_τ′ = Int(round(1 / (T′ * Δτ)))
    global n_τ = n_τ′

    global σ = ones_Z2_gauge_field_DPI(Int, lattice, n_τ)
    global model_Z2 = IsingGaugeTheoryDPIMetropolisMC(Float64, σ, g, h_Z2, Δτ)
    global model_Z2_fermion = Z2SpinlessFermionSimpleDQMC(Float64, σ, n_wrap, t, Δτ)
    global aux = Z2SpinlessFermionSimpleAuxField(Float64, model_Z2_fermion, σ)

    global s = ones_Ising_field_DPI(Int, lattice, n_τ)
    global model_TFIM = TransverseFieldIsingModelDPIMetropolisMC(Float64, s, 0.0, h_TFIM, Δτ)

    global model_coupling = Z2IsingMinimalCouplingDPI{Float64}(J, Δτ)
end

###################################################################################
# Heating up.
###################################################################################

function sweep()
    for τ in 1 : n_τ
        for i in lattice_sites
            pure_TFIM_weight_ratio = weight_ratio(model_TFIM, s, i, τ)
            coupling_weight_ratio = weight_ratio_s(model_coupling, σ, s, i, τ)
            accept_rate = pure_TFIM_weight_ratio * coupling_weight_ratio

            if rand() < accept_rate
                update!(model_TFIM, s, i, τ)
            end
        end        
    end

    for τ in 1 : n_τ - 1
        for b in lattice_bonds
            pure_Z2_weight_ratio = weight_ratio(model_Z2, σ, b, τ)
            Z2_Ising_coupling_weight_ratio = weight_ratio_σ(model_coupling, σ, s, b, τ)
            Z2_fermion_weight_ratio = weight_ratio(model_Z2_fermion, aux, b, τ)
            accept_rate = pure_Z2_weight_ratio * Z2_Ising_coupling_weight_ratio * Z2_fermion_weight_ratio

            if rand() < accept_rate
                update!(model_Z2_fermion, aux, b, τ)
            end
        end
        propagate_forward!(model_Z2_fermion, aux, τ)
    end

    for τ in n_τ : -1 : 2
        for b in lattice_bonds
            pure_Z2_weight_ratio = weight_ratio(model_Z2, σ, b, τ)
            coupling_weight_ratio = weight_ratio_σ(model_coupling, σ, s, b, τ)
            accept_rate = pure_Z2_weight_ratio * coupling_weight_ratio

            if rand() < accept_rate
                update!(model_Z2_fermion, aux, b, τ)
            end
        end
        propagate_backward!(model_Z2_fermion, aux, τ)
    end
end

function heatup()
    println("Heating up.")

    progress = Progress(n_heat)

    for _ in n_heat
        sweep()
        if show_progress
            next!(progress)
        end
    end
end

###################################################################################
# Observe.
###################################################################################

function observe()
    println("Starting observation.")

    flux_history = []
    magnetization_history = []

    progress = Progress(n_bin * n_sweep)

    for _ in 1 : n_bin
        flux_bin = []
        magnetization_bin = []
        for _ in 1 : n_sweep
            sweep()
            push!(flux_bin, flux_average(σ, n_τ))
            push!(magnetization_bin, magnetization(s, n_τ))
            if show_progress
                next!(progress)
            end
        end
        push!(flux_history, mean(flux_bin))
        push!(magnetization_history, mean(magnetization_bin))
    end

    println("Calculation completed.")

    return mean(flux_history), std(flux_history), mean(magnetization_history), std(magnetization_history)
end

###################################################################################
# Sweeping throught the phase diagram.
###################################################################################

digits_kept = 10
standard_output(data) = rpad(string(round(data, digits = digits_kept)), digits_kept + 2, "0")

T_range = LinRange(0.5, 3.0, 50)
h_Z2_range = LinRange(0.1, 2, 40)

open(working_path * "flux", "w") do file
    println(file, "Average Z2 flux and standard variation")
    println(file, "T             h             flux          std")
    println(file, "========================================================")
end

open(working_path * "magnetization", "w") do file
    println(file, "Average Ising magnetization and standard variation")
    println(file, "T             h             mag           std")
    println(file, "========================================================")
end

for T′ in T_range
    for h_Z2′ in h_Z2_range
        refresh_params(t, h_Z2′, h_TFIM, J, g, T′)
        heatup()
        flux_mean, flux_std, mag_mean, mag_std = observe()
         
        output_flux = "$(standard_output(T′))  $(standard_output(h_Z2′))  "
        output_flux *= "$(standard_output(flux_mean))  $(standard_output(flux_std))"
        output_mag = "$(standard_output(T′))  $(standard_output(h_Z2′))  "
        output_mag *= "$(standard_output(mag_mean))  $(standard_output(mag_std))"
        
        open(working_path * "flux", "a") do file
            println(file, output_flux)
        end

        open(working_path * "magnetization", "a") do file
            println(file, output_mag)
        end
    end
end