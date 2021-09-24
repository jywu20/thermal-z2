struct SpinlessFermionDQMCConfiguration2D{F <: Float64, I <: Integer, L <: Integer}
    n_side::L
    G::Matrix{F}
    s::Array{I}    
end

struct SpinlessFermion 
    
end

struct SquareLattice2D
end