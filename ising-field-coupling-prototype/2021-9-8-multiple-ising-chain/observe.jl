function magnetization(;kx, ky, nx, ny, n_bin, n_heat, n_sweep)
    sim_params_heat = AnisotropicIsing2DSimParams(nx, ny, n_heat)
    ham_params = AnisotropicIsing2DHamParams(kx, ky)
    ising_field = init(sim_params_heat)
    run!(ham_params, sim_params_heat, ising_field)

    sim_params = AnisotropicIsing2DSimParams(nx, ny, n_sweep)
    mag, _ = binning(() -> run!(ham_params, sim_params, ising_field; observe = abs ∘ mean), n_bin)

    mag
end