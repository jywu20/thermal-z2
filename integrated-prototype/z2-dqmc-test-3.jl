using Plots
using Statistics
using ProgressMeter
using LaTeXStrings

##

include("lib.jl")

##

n_sweep = 5
n_side = 4
n_τ = 10
n_wrap = 10
t = 1.0
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

G_err_history = []
weight_ratio_err_history = []

progress = Progress(n_sweep * 2 * (n_τ - 1) * 2 * n_side^2)

for _ in 1 : n_sweep
    for τ in 1 : n_τ - 1
        for b in 1 : 2 * n_side^2
            update!(model, aux, b, τ)
            push!(G_err_history, relative_err(G_τ_τ_def(model, σ, τ), aux.G[:, :, τ]))
            push!(weight_ratio_err_history, relative_err(weight_ratio(model, aux, b, τ), weight_ratio_def(model, σ, b, τ)))
            next!(progress)
        end
        propagate_forward!(model, aux, τ)
    end
    
    for τ in n_τ : -1 : 2
        for b in 1 : 2 * n_side^2
            update!(model, aux, b, τ)
            push!(G_err_history, relative_err(G_τ_τ_def(model, σ, τ), aux.G[:, :, τ]))
            push!(weight_ratio_err_history, relative_err(weight_ratio(model, aux, b, τ), weight_ratio_def(model, σ, b, τ)))
            next!(progress)
        end
        propagate_backward!(model, aux, τ)
    end
end

plot(G_err_history, legend = false)

##

plot(weight_ratio_err_history, legend = false)

##

τ = 10
relative_err(G_τ_τ(model, aux, τ), G_τ_τ_def(model, σ, τ))

##

τ = 10
relative_err(G_τ_τ(model, aux, τ), aux.G[:, :, τ])
G_τ_τ!(model, aux, τ)
relative_err(G_τ_τ_def(model, σ, τ), aux.G[:, :, τ])