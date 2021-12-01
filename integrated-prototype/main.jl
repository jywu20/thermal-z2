using Plots
using Statistics
using ProgressMeter
using LaTeXStrings

include("lib.jl")

n_sweep = 100
n_heat = 100
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

lattice = SquareLattice2DPeriodic(n_side)
lattice_sites = sites(lattice)
lattice_bonds = bonds(lattice)

σ = ones_Z2_gauge_field_DPI(Int, lattice, n_τ)
model_Z2 = IsingGaugeTheoryDPIMetropolisMC(Float64, σ, g, h_Z2, Δτ)
aux = Z2SpinlessFermionSimpleAuxField(Float64, model_Z2, σ)

s = ones_Ising_field_DPI(Int, lattice, n_τ)
model_TFIM = TransverseFieldIsingModelDPIMetropolisMC(Float64, s, 0.0, h_TFIM, Δτ)

for _ in n_heat
    for τ in 1 : n_τ
        for i in lattice_sites
            pure_TFIM_weight_ratio = weight_ratio(model_TFIM, s, i, τ)
            # TODO: coupling 
        end        
    end

    for τ in 1 : n_τ
        
    end
end