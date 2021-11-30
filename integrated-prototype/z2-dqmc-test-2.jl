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
t = 8.0
Δτ = 0.01
β = n_τ * Δτ 

lattice = SquareLattice2DPeriodic(n_side)
σ = ones_Z2_gauge_field_DPI(Int, lattice, n_τ)
σ[12, 1] = σ[14, 1] = -1
σ[2, 2] = σ[3, 2] = σ[1, 2] = σ[7, 2] = -1
σ[11, 3] = σ[6, 3] = σ[18, 3] = σ[24, 3] = -1
σ[10, 4] = σ[1, 4] = -1
σ[8, 5] = -1
σ[1, 6] = σ[2, 6] = σ[3, 6] = σ[15, 6] = σ[14, 6] = -1

model = Z2SpinlessFermionSimpleDQMC(Float64, σ, n_wrap, t, Δτ)
aux = Z2SpinlessFermionSimpleAuxField(Float64, model, σ)

##

τ = 10

relative_err(aux.G[:, :, τ], G_τ_τ_def(model, σ, τ))

##

τ = 2
V, D, U = B_β_τ_vdu(model, aux, τ)
relative_err(B_β_τ_def(model, σ, τ), V * D * U)

##

τ = 2
U, D, V = B_τ_0_udv(model, aux, τ)
relative_err(B_τ_0_def(model, σ, τ), U * D * V)

##

b = 10
τ = 2 

relative_err(weight_ratio(model, aux, b, τ), accept_rate_def(model, aux.σ, b, τ))

##

heatmap(I + B_β_τ_def(model, σ, 0))

##
σ′ = copy(σ)
σ′[b, τ] *= -1
heatmap(I + B_β_τ_def(model, σ′, 0))

##

b = 12
τ = 6
update!(model, aux, b, τ)
(relative_err(G_τ_τ_def(model, aux.σ, τ), aux.G[:, :, τ]), 
relative_err(weight_ratio(model, aux, b, τ), accept_rate_def(model, aux.σ, b, τ)))

##
σ′ = copy(σ)
B = B_τ_0_def(model, σ, τ)
σ′[b, τ] *= -1
B′ = B_τ_0_def(model, σ′, τ)
@show relative_err(B′, (I + Δ_mat(model, aux, b, τ)) * B)
@show relative_err(B, (I + Δ_mat_def(model, σ′, b, τ)) * (I + Δ_mat(model, aux, b, τ)) * B)

##

τ = 1
G_old = aux.G[:, :, τ + 1]
propagate_forward!(model, aux, τ)
G_new = aux.G[:, :, τ + 1]
relative_err(G_old, G_new)

##

τ = 10 
G_old = aux.G[:, :, τ - 1]
propagate_backward!(model, aux, τ)
G_new = aux.G[:, :, τ - 1]
relative_err(G_old, G_new)