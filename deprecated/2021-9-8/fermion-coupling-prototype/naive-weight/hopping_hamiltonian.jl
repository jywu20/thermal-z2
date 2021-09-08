using Plots

##

function back_into_range(idx, upper)
    if idx > upper
        return idx % upper
    end
    (idx - upper) % upper + upper
end

struct Z2Gauge2DHamParams
    t::Float64
end

struct Z2GaugeField2D 
    sub1::Matrix{Int64}
    sub2::Matrix{Int64}
end

const Site2D = Tuple{Int64, Int64}

function getindex(z2_field::Z2GaugeField2D, i::Site2D, j::Site2D)::Float64
    ix, iy = i
    jx, jy = j
    
    if ix != jx && iy != jy
        return 0.0
    end

    if ix == jx
        if iy == jy + 1
            return z2_field.sub1[ix, jy]
        end
        if jy == iy + 1
            return z2_field.sub1[ix, iy]
        end
    end
    
    if iy == jy
        if ix == jx + 1
            return z2_field.sub2[jx, iy]
        end
        if jx == ix + 1
            return z2_field.sub2[ix, iy]
        end
    end
    return 0.0
end

function getindex(z2_field::Z2GaugeField2D, i::Site2D, j::Site2D)
    
end

struct Z2Gauge2DSimParams 
    n_side::Int64
end

##

function hopping_ham(ham_params::Z2Gauge2DHamParams, sim_params::Z2Gauge2DSimParams)
    n_side = sim_params.n_side
    t = ham_params.t

    T = zeros((n_side, n_side))


end