using Plots
using JLD

##

magnetization, magnetic_susceptibility, T_range, h_range = jldopen("./data.jld", "r") do file
    (read(file, name) for name in ["magnetization", "magnetic_susceptibility", "T_range", "h_range"])
end

##

p = heatmap(h_range, T_range, magnetization, xlabel = L"h", ylabel = L"T")
ylims!(p, T_range[1], T_range[end])
xlims!(p, h_range[1], h_range[end])
 
##

p = heatmap(h_range[1 : end - 1], T_range, magnetic_susceptibility, xlabel = L"h", ylabel = L"T")
ylims!(p, T_range[1], T_range[end])
xlims!(p, h_range[1], h_range[end - 1])