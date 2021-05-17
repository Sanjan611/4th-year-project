using RHEOS
using Plots

include("springpot_helper.jl")

# βs = [0.1, 0.3, 0.5, 0.7, 0.9]

# for (i, β) in enumerate(βs)
#     # Model
#     model_params = (cᵦ = 1.8, β = β)
#     model = RheoModel(Springpot, model_params)
#     model_i = RheoModel(Springpot_i, model_params)

#     time_sim = timeline(t_start = 0, t_end = 8, step = 0.02)
#     load_sim = strainfunction(time_sim, t->f(t,model_params.β));

#     σ = σA(model, load_sim)

#     t = time_sim.t

#     e1 = abs.(σ .- σB(model, load_sim))
#     e2 = abs.(σ .- σR(model_i, load_sim))

#     # Plotting and saving image
#     plt = plot(t, [e1, e2], label=["B" "R"], linewidth=2.0, legendfontsize = 10)
#     xlabel!("time point")
#     ylabel!("absolute error")
#     title!(string(β))
#     savefig("images/comp_deriv_beta_" * string(10*β) * ".svg")
#     display(plt)

# end

## NOW I WANT TO COMPARE ERROR OVER A LOT OF BETAS AT SINGLE TIME POINT 

βs = range(0.01, 0.90, length = 30)
errors = [[], []]

time_sim = timeline(t_start = 0, t_end = 8, step = 0.02)

for (i, β) in enumerate(collect(βs))

    # Model
    model_params = (cᵦ = 1.8, β = β)
    model = RheoModel(Springpot, model_params)
    model_i = RheoModel(Springpot_i, model_params)

    load_sim = strainfunction(time_sim, t->f(t,model_params.β));
    σ = σA(model, load_sim)

    # Errors
    e1 = abs.(σ .- σB(model, load_sim)) # error with original boltzmann
    e2 = abs.(σ .- σR(model_i, load_sim)) # error with ramp boltzmann

    # consider only the last time point
    push!(errors[1], e1[5])
    push!(errors[2], e2[5])

end

t = time_sim.t
plt = plot(collect(βs), errors, label = ["B" "R"], linewidth = 2.0, legendfontsize=10)
xlabel!("β")
ylabel!("error at last time point")
title!("Error at last time point vs β")
savefig("images/comp_deriv_beta_one_time.svg")
display(plt)
