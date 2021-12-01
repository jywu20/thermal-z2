using Statistics

include("lib.jl")

###################################################################################
# Parameters
###################################################################################

n_sweep = 100
n_heat = 100
n_bin = 10
n_side = 4
n_τ = 10
n_wrap = 10

Δτ = 0.01
β = n_τ * Δτ 

t = 1.0
h_Z2 = 1.0
h_TFIM = 1.0
J = 1.0
g = 1.0

###################################################################################
# Lattice and fields initialization.
###################################################################################

lattice = SquareLattice2DPeriodic(n_side)
lattice_sites = sites(lattice)
lattice_bonds = bonds(lattice)

σ = ones_Z2_gauge_field_DPI(Int, lattice, n_τ)
model_Z2 = IsingGaugeTheoryDPIMetropolisMC(Float64, σ, g, h_Z2, Δτ)
aux = Z2SpinlessFermionSimpleAuxField(Float64, model_Z2, σ)

s = ones_Ising_field_DPI(Int, lattice, n_τ)
model_TFIM = TransverseFieldIsingModelDPIMetropolisMC(Float64, s, 0.0, h_TFIM, Δτ)

model_coupling = Z2IsingMinimalCouplingDPI{Float64}(J, Δτ)

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
            coupling_weight_ratio = weight_ratio_σ(model_coupling, σ, s, i, τ)
            accept_rate = pure_Z2_weight_ratio * coupling_weight_ratio

            if rand() < accept_rate
                update!(model_Z2, aux, b, τ)
            end
        end
        propagate_forward!(model_Z2, aux, τ)
    end

    for τ in n_τ : -1 : 2
        for b in lattice_bonds
            pure_Z2_weight_ratio = weight_ratio(model_Z2, σ, b, τ)
            coupling_weight_ratio = weight_ratio_σ(model_coupling, σ, s, i, τ)
            accept_rate = pure_Z2_weight_ratio * coupling_weight_ratio

            if rand() < accept_rate
                update!(model_Z2, aux, b, τ)
            end
        end
        propagate_backward!(model_Z2, aux, τ)
    end
end

for _ in n_heat
    sweep()
end

###################################################################################
# Observe.
###################################################################################

flux_average_history = []
for _ in 1 : n_bin
    flux_bin = []
    for _ in 1 : n_sweep
        sweep()
        push!(flux_bin, flux_average(σ, n_τ))
    end
    push!(flux_average_history, mean(flux_bin))
end