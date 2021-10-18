function back_into_range(idx, upper)
    if idx > upper
        return idx % upper
    end
    (idx - upper) % upper + upper
end

commutator(a, b) = a * b - b * a

relative_err(a, b) = norm(a - b) / norm(a)