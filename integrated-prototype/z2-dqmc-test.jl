using Plots
using Statistics
using ProgressMeter
using LaTeXStrings

##

include("lib.jl")

##

n_side = 4
n_τ = 10
lattice = SquareLattice2DPeriodic(n_side)
σ = ones_Z2_gauge_field_DPI(Int, lattice, n_τ)
σ[12, 1] = -1
σ[14, 1] = -1
σ[2, 2] = -1

model = Z2SpinlessFermionDQMC{Float64, Matrix{Float64}}(1.0, 0.5, 4)

τ = 1

##
@show T_def(model, σ, τ)[12, 16]
@show T_def(model, σ, τ)[1, 5]
@show T_def(model, σ, τ)[2, 14]
@show T_def(model, σ, τ)[14, 15]

# Nothing wrong is found.

##

Δτ_list = collect(0.001:0.002:0.05)
progress = Progress(length(Δτ_list))
error_1_list = []
error_2_list = []

for Δτ in Δτ_list
    τ = 1
    model = Z2SpinlessFermionDQMC{Float64, Matrix{Float64}}(1.0, Δτ, 4)

    n_side = 4
    n_site = n_side^2
    n_τ = 10
    lattice = SquareLattice2DPeriodic(n_side)
    σ = random_Z2_gauge_field_DPI(Int, lattice, n_τ)

    B1 = B_def(model, σ, τ)
    B2 = map(b -> B_at_bond(model, σ, b, τ), bonds(lattice)) |> prod
    push!(error_1_list, relative_err(B1, B2))
    push!(error_2_list, relative_err(B1 - I, B2 - I))

    next!(progress)
end

p = plot(Δτ_list, error_1_list, legend = false)
plot!(p, Δτ_list, error_2_list, legend = false)
xlabel!(p, L"\Delta\tau")
ylabel!(p, "relative error")

##
