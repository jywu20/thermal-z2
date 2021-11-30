# Calculate B matrices.

B_mat(model::Z2SpinlessFermionSimpleDQMC, aux::Z2SpinlessFermionSimpleAuxField, τ) = B_def(model, aux.σ, τ)

function B_τ_0(model::Z2SpinlessFermionSimpleDQMC, aux::Z2SpinlessFermionSimpleAuxField, τ)
    σ = aux.σ
    n_site = site_number(σ.lattice)
    B = Matrix(I, n_site, n_site)
    for τ′ in 1 : τ
        B = B_mat(model, σ, τ′) * B        
    end
    B
end

function B_β_τ(model::Z2SpinlessFermionSimpleDQMC, aux::Z2SpinlessFermionSimpleAuxField, τ)
    σ = aux.σ
    n_site = site_number(σ.lattice)
    n_τ = σ.n_τ
    B = Matrix(I, n_site, n_site)
    for τ′ in n_τ : τ + 1
        B = B * B_def(model, σ, τ′)
    end
    B
end