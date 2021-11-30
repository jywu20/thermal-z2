using Plots
using Statistics
using ProgressMeter
using LaTeXStrings

##

include("lib.jl")

##

n_side = 4
n_τ = 10
n_wrap = 10
t = 1.0
Δτ = 0.5
β = n_τ * Δτ 

lattice = SquareLattice2DPeriodic(n_side)
σ = ones_Z2_gauge_field_DPI(Int, lattice, n_τ)
σ[12, 1] = -1
σ[14, 1] = -1
σ[2, 2] = -1

model = Z2SpinlessFermionSimpleDQMC(Float64, σ, n_wrap, t, Δτ)
aux = Z2SpinlessFermionSimpleAuxField(Float64, model, σ)

##

τ = 2

relative_err(aux.G[:, :, 1], G_τ_τ_def(model, σ, 1))
# Error here