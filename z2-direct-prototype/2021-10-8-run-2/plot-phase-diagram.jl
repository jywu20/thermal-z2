using Plots
using LaTeXStrings

##

T_steps = 80
T_range = LinRange(0.1, 4.1, T_steps)

h_steps = 80
h_range = LinRange(0.1, 3.1, h_steps)
Δh = (h_range[end] - h_range[1]) / h_steps

# Since the computation is not finished, we have to tailor `T_range`
T_range = T_range[1:end-1]

##

# Since we are going to run the program by `julia main.jl`, the working directory can be set to "./"
working_path = "D:\\Projects\\thermal-z2\\z2-direct-prototype\\2021-10-8-run-2\\" 
output_text_name = "output"
mag = open(working_path * output_text_name) do file
    mag_str = readlines(file)[8:end]
    mag = map(x -> parse(Float64, split(x)[end]), mag_str)
    # Since the computation is not finished, we have to tailor the output data
    mag = mag[1:length(h_range) * length(T_range)]
    reshape(mag, length(h_range), length(T_range))'
end

##

heatmap(h_range, T_range, mag, xlabel = L"h", ylabel = L"T")

##

magnetic_susceptibility = [abs(mag[i, j + 1] - mag[i, j]) / Δh for i in 1 : length(T_range), j in 1 : length(h_range) - 1]
p = heatmap(h_range[1 : end - 1], T_range, magnetic_susceptibility, xlabel = L"h", ylabel = L"T")
ylims!(p, T_range[1], T_range[end])
xlims!(p, h_range[1], h_range[end - 1])