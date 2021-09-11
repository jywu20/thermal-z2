using Plots
using LaTeXStrings

J_range = LinRange(0.0, 4.0, 40)
Δτ_range = LinRange(0.001, 0.1, 20)

heatmap(J_range, Δτ_range, [Δτ * J for Δτ in Δτ_range, J in J_range])

##

h = 2.0
heatmap(J_range, Δτ_range, [atanh(Δτ * h) for Δτ in Δτ_range, J in J_range])

##

h = 2.0
heatmap(J_range, Δτ_range, [atanh(Δτ * h) / (Δτ * J) for Δτ in Δτ_range, J in J_range])

##

h = 2.0
heatmap(J_range, Δτ_range, [(Δτ * J) / atanh(Δτ * h) for Δτ in Δτ_range, J in J_range], xlabel = L"J", ylabel = L"\Delta \tau", title = L"J_{xy} / J_\tau")

##

function heatmap_jxy_jtau_ratio(h)
    fig_title = "J_{xy} / J_\\tau \\ \\mathrm{with} \\ h = $h"
    heatmap(J_range, Δτ_range, [(Δτ * J) / atanh(Δτ * h) for Δτ in Δτ_range, J in J_range], xlabel = L"J", ylabel = L"\Delta \tau", title = latexstring(fig_title))
end

function contourf_jxy_jtau_ratio(h)
    fig_title = "J_{xy} / J_\\tau \\ \\mathrm{with} \\ h = $h"
    contourf(J_range, Δτ_range, [(Δτ * J) / atanh(Δτ * h) for Δτ in Δτ_range, J in J_range], xlabel = L"J", ylabel = L"\Delta \tau", title = latexstring(fig_title))
end

##
contourf_jxy_jtau_ratio(5.0)
