using ProgressMeter
include("z2.jl")

##

n_heat = 100
n_steps = 20
n_side = 10

J = 1.0
h = 1.0
T = 1.0
β = 1 / T

lattice = PeriodicSquareLattice2DWithBonds(n_side)
z2field = all_ones_discrete_path_integral_z2_gauge_field(lattice, n_steps)
params = Z2GaugeTheoryDiscretePathIntegralMetropolisParams(J, h, n_side, n_steps, β)

progress = Progress(n_heat)

sweep!(z2field, params, n_heat; observe = _ -> next!(progress))