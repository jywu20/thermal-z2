using Plots
using LaTeXStrings

J_range = LinRange(0.2, 4.0, 100)
Δτ_range = LinRange(0.01, 1.0, 100)

heatmap(J_range, Δτ_range, [Δτ * J for Δτ in Δτ_range, J in J_range])

##

h = 2.0
heatmap(J_range, Δτ_range, [atanh(exp(- 2 * Δτ * h)) for Δτ in Δτ_range, J in J_range])

##

h = 2.0
heatmap(J_range, Δτ_range, [atanh(exp(- 2 * Δτ * h)) / (Δτ * J) for Δτ in Δτ_range, J in J_range])

##

h = 2.0
heatmap(J_range, Δτ_range, [(Δτ * J) / atanh(exp(- 2 * Δτ * h)) for Δτ in Δτ_range, J in J_range], xlabel = L"J", ylabel = L"\Delta \tau", title = L"J_{xy} / J_\tau")

##

function heatmap_jxy_jtau_ratio(h)
    fig_title = "J_{xy} / J_\\tau \\ \\mathrm{with} \\ h = $h"
    heatmap(J_range, Δτ_range, [(Δτ * J) / atanh(exp(- 2 * Δτ * h)) for Δτ in Δτ_range, J in J_range], xlabel = L"J", ylabel = L"\Delta \tau", title = latexstring(fig_title))
end

function contourf_jxy_jtau_ratio(h)
    fig_title = "J_{xy} / J_\\tau \\ \\mathrm{with} \\ h = $h"
    contourf(J_range, Δτ_range, [(Δτ * J) / atanh(exp(- 2 * Δτ * h)) for Δτ in Δτ_range, J in J_range], xlabel = L"J", ylabel = L"\Delta \tau", title = latexstring(fig_title))
end

heatmap_jxy_jtau_ratio(1.0)
#contourf_jxy_jtau_ratio(1.0)

##

function heatmap_jxy_jtau_ratio_log(h)
    fig_title = "\\log(J_{xy} / J_\\tau) \\ \\mathrm{with} \\ h = $h"
    p = heatmap(J_range, Δτ_range, [log((Δτ * J) / atanh(exp(- 2 * Δτ * h))) for Δτ in Δτ_range, J in J_range], xlabel = L"J", ylabel = L"\Delta \tau", title = latexstring(fig_title))
    p
end

function contourf_jxy_jtau_ratio_log(h)
    fig_title = "\\log(J_{xy} / J_\\tau) \\ \\mathrm{with} \\ h = $h"
    p = contourf(J_range, Δτ_range, [log((Δτ * J) / atanh(exp(- 2 * Δτ * h))) for Δτ in Δτ_range, J in J_range], xlabel = L"J", ylabel = L"\Delta \tau", title = latexstring(fig_title))
    p
end

contourf_jxy_jtau_ratio_log(0.001)

##

function contourf_jxy_jtau_2_ratio_log(h)
    fig_title = "\\log(J_{xy} / J_\\tau) \\ \\mathrm{with} \\ h = $h"
    p = contourf(J_range, Δτ_range, [log((Δτ * J) / (- log(tanh(Δτ * h)) / 2)) for Δτ in Δτ_range, J in J_range], xlabel = L"J", ylabel = L"\Delta \tau", title = latexstring(fig_title))
    p
end

contourf_jxy_jtau_2_ratio_log(10)