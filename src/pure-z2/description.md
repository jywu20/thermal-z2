# Description of files in this folder

- `anisotropic_ising.jl` is the core of simulations. 
  It contains a cluster-update algorithm for 2D anisotropic Ising model.
- `observe.jl` contains observables for anisotropic Ising model.
- `utils.jl` contains utilities, including a function imposing the periodic boundary condition 
  and a function for binning.
- `magnetic-ordering.jl` investigates magnetic orders of a 1D transverse field Ising chain, which is 
  mapped into a 2D anosotropic Ising model. The results agree with the theoretical derived conclusion that 
  there is no magnetic ordering in 1D transverse field Ising chain under finite temperature.
  