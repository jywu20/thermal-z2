abstract type AbstractModel end

abstract type AbstractConfigration end

"""
Abstract function for doing one Metropolis step.
"""
function metropolis!(::M, ::C) where {M <: AbstractModel, C <: AbstractConfigration}
    error("Metropolis updating not defined for model $M on configuration type $C.")
end

"""
Abstract function for doing one Wolff step.
"""
function wolff!(::M, ::C) where {M <: AbstractModel, C <: AbstractConfigration}
    error("Wolff cluster updating not defined for model $M on configuration type $C.")
end