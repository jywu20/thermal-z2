function B_at_bond(model::Z2SpinlessFermionDQMC{F, M}, σ::Z2GaugeFieldDPI, b, τ) where {
    F <: AbstractFloat, M <: AbstractMatrix
}
    lattice = σ.lattice
    n_site = site_number(lattice)
    i, j = bond_to_sites(lattice, b)
    t = model.t
    Δτ = model.Δτ

    T = zeros(F, n_site, n_site)
    T[i, j] = T[j, i] = Δτ * t * σ[b, τ]
    M(exp(T))
end

B_at_bond(model::Z2SpinlessFermionDQMC, config::Z2SpinlessFermionAuxField, b, τ) = B_at_bond(model, config.σ, b, τ)

function B_prod_bonds(model::Z2SpinlessFermionDQMC{F, M}, σ::Z2GaugeFieldDPI, τ) where {
    F <: AbstractFloat, M <: AbstractMatrix
}
    lattice = σ.lattice
    map(b -> B_at_bond(model, σ, b, τ), bonds(lattice)) |> prod
end