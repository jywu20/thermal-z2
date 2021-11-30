# Definitions used in DQMC, with poor performance, used for benchmark

"""
B-matrices by definition, without any numerical acceleration.
"""
function B_def(model::Z2SpinlessFermionSimpleDQMC, σ::Z2GaugeFieldDPI, τ)
    Δτ = model.Δτ
    t = model.t
    n_site = site_number(σ.lattice)
    T = zeros(n_site, n_site)
    for b in bonds(lattice)
        i, j = bond_to_sites(lattice, b)
        T[i, j] = T[j, i] = σ[b, τ]
    end
    exp(Δτ * t * T)
end

function B_τ_0_def(model::Z2SpinlessFermionSimpleDQMC, σ::Z2GaugeFieldDPI, τ)
    n_site = site_number(σ.lattice)
    B = Matrix(I, n_site, n_site)
    for τ′ in 1 : τ
        B = B_def(model, σ, τ′) * B        
    end
    B
end

function B_β_τ_def(model::Z2SpinlessFermionSimpleDQMC, σ::Z2GaugeFieldDPI, τ)
    n_site = site_number(σ.lattice)
    n_τ = σ.n_τ
    B = Matrix(I, n_site, n_site)
    for τ′ in n_τ : τ + 1
        B = B * B_def(model, σ, τ′)
    end
    B
end

function G_τ_τ_def(model::Z2SpinlessFermionSimpleDQMC, σ::Z2GaugeFieldDPI, τ)
    inv(I + B_τ_0_def(model, σ, τ) * B_β_τ_def(model, σ, τ))
end

function T_def(model::Z2SpinlessFermionSimpleDQMC, σ::Z2GaugeFieldDPI, τ)
    t = model.t
    n_site = site_number(σ.lattice)
    T = zeros(n_site, n_site)
    for b in bonds(lattice)
        i, j = bond_to_sites(lattice, b)
        T[i, j] = T[j, i] = σ[b, τ]
    end
    - t * T
end