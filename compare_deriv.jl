using RHEOS
using Plots

include("springpot_helper.jl")

# Model
model_params = (cᵦ = 1.8, β = 0.5)
model = RheoModel(Springpot, model_params)
model_i = RheoModel(Springpot_i, model_params)

# time and strain data
time_sim = timeline(t_start = 0, t_end = 10, step = 0.01)
load_sim = strainfunction(time_sim, t->f(t,model_params.β));

# Get the absolute error using the 4 different methods to calculate the derivatives
σ = σA(model, load_sim)
t = time_sim.t
e1 = abs.(σ .- σB(model, load_sim))
e2 = abs.(σ .- σR(model_i, load_sim))
e3 = abs.(σ .- σL1(model, load_sim))
e4 = abs.(σ .- σL12(model, load_sim))

# Plot
plt = plot(t, [e1, e2, e3, e4], label=["B" "R" "L1" "L12"], 
        linewidth=2, legendfontsize = 10)
xlabel!("time point")
ylabel!("absolute error")
title!("Error for all methods over time")
savefig("images/abs_error_4_methods.svg")
display(plt)

