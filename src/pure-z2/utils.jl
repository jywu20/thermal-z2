function back_into_range(idx, upper)
    if idx > upper
        return idx % upper
    end
    (idx - upper) % upper + upper
end

function binning(func, n_bin; verbose = false)
    if verbose
        println("$n_bin bins are dealt with.")
    end
    bins = []
    for bin_count in 1 : n_bin
        this_bin_result = mean(func())
        if verbose
            println("Bin $bin_count is finished.")
        end
        push!(bins, this_bin_result)
    end

    (mean(bins), std(bins))
end