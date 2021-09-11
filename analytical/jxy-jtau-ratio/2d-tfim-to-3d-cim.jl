using Plots
using LaTeXStrings

function heatmap_jxy_jtau_ratio(h, J_range, Δτ_range; notitle = false)
    fig_title = "J_{xy} / J_\\tau \\ \\mathrm{with} \\ h = $h"
    if notitle
        return heatmap(J_range, Δτ_range, [(Δτ * J) / atanh(Δτ * h) for Δτ in Δτ_range, J in J_range], xlabel = L"J", ylabel = L"\Delta \tau")
    end
    heatmap(J_range, Δτ_range, [(Δτ * J) / atanh(Δτ * h) for Δτ in Δτ_range, J in J_range], xlabel = L"J", ylabel = L"\Delta \tau", title = latexstring(fig_title))
end

function contourf_jxy_jtau_ratio(h, J_range, Δτ_range; notitle = false)
    fig_title = "J_{xy} / J_\\tau \\ \\mathrm{with} \\ h = $h"
    if notitle
        return contourf(J_range, Δτ_range, [(Δτ * J) / atanh(Δτ * h) for Δτ in Δτ_range, J in J_range], xlabel = L"J", ylabel = L"\Delta \tau")
    end
    contourf(J_range, Δτ_range, [(Δτ * J) / atanh(Δτ * h) for Δτ in Δτ_range, J in J_range], xlabel = L"J", ylabel = L"\Delta \tau", title = latexstring(fig_title))
end

##

J_range = LinRange(0.0, 4.0, 40)
Δτ_range = LinRange(0.001, 0.1, 20)
h = 2.0

contourf_jxy_jtau_ratio(h, J_range, Δτ_range, notitle = true)
