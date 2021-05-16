using CSV
using DataFrames
using Plots

df = CSV.read("springpot_1_withGaussianNoise.csv", DataFrame)
labels = ["A", "B", "R", "L1", "L12"]
error = df.error
error_m = reshape(error, (5,5))
timetaken = df.time_total
timetaken_m = reshape(timetaken, (5,5))

# function f(x, y)
#     r = sqrt(x^2 + y^2)
#     return cos(r) / (1 + r)
# end
# x = range(0, 2π, length = 30)
# heatmap(x, x, f, c = :autumn1)

# heatmap(labels, labels, error_m, c = :haline)]
# title!("Error")
# xlabel!("σ")
# ylabel!("How variable part is calculated")

heatmap(labels, labels, timetaken_m, c = :haline)
title!("Time taken")
xlabel!("σ")
ylabel!("How variable part is calculated")