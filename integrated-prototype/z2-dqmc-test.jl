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
# Trotter decomposition test: whether the approximation of B_τ by ∑_{⟨i,j⟩} B_{ijτ} is acceptable.

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
    B2 = B_prod_bonds(model, σ, τ)
    push!(error_1_list, relative_err(B1, B2))
    push!(error_2_list, relative_err(B1 - I, B2 - I))

    next!(progress)
end

p = plot(Δτ_list, error_1_list, legend = false)
plot!(p, Δτ_list, error_2_list, legend = false)
xlabel!(p, L"\Delta\tau")
ylabel!(p, "relative error")

##
# Trotter decomposition test: whether updating B matrices via Δ matrices is acceptable.

Δτ_list = collect(0.001:0.002:0.05)
error_list = []

for Δτ in Δτ_list
    model = Z2SpinlessFermionDQMC{Float64, Matrix{Float64}}(1.0, Δτ, 4)

    n_side = 4
    n_site = n_side^2
    n_τ = 10
    lattice = SquareLattice2DPeriodic(n_side)
    σ = random_Z2_gauge_field_DPI(Int, lattice, n_τ)

    τ = 4
    b = 10
    B_τ_this = B_def(model, σ, τ)
    σ′ = copy(σ)
    σ′[b, τ] *= -1
    B′_τ_this = B_def(model, σ′, τ)

    push!(error_list, relative_err((I + Δ_mat(model, σ, b, τ)) * B_τ_this, B′_τ_this))
end

plot(Δτ_list, error_list, legend = false, xlabel=L"\Delta\tau")

##
# Test the accumulation of Trotter error.

Δτ = 0.005
model = Z2SpinlessFermionDQMC{Float64, Matrix{Float64}}(1.0, Δτ, 4)

n_side = 4
n_site = n_side^2
n_τ = 10
lattice = SquareLattice2DPeriodic(n_side)
σ = random_Z2_gauge_field_DPI(Int, lattice, n_τ)

τ = 4

b_list = bonds(lattice)
error_list = []

B_τ_this = B_def(model, σ, τ)

for b in b_list
    σ′ = copy(σ)
    σ′[b, τ] *= -1
    B′_τ_this = B_def(model, σ′, τ)

    copy!(B_τ_this, (I + Δ_mat(model, σ, b, τ)) * B_τ_this)

    push!(error_list, relative_err(B_τ_this, B′_τ_this))
end

plot(b_list, error_list, legend = false)

##
