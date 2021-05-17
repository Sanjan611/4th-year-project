using RHEOS
using Plots

include("springpot_helper.jl")

# Model
model_params = (cᵦ = 1.8, β = 0.3)
model = RheoModel(Springpot, model_params)
model_i = RheoModel(Springpot_i, model_params)

type = ["B", "R", "L1", "L12"]
σs = [σB, σR, σL1, σL12]

stepsizes = [0.01, 0.03, 0.05, 0.07, 0.09]

for (i, σi) in enumerate(σs)
    # Iterate over the 4 different derivative calc methods

    errors = []
    timepoints = []

    for (j,ss) in enumerate(stepsizes)
        # Iterate over the step sizes

        time_sim = timeline(t_start = 0, t_end = 10, step = ss)
        load_sim = strainfunction(time_sim, t->f(t,model_params.β));

        σ = σA(model, load_sim)
        t = time_sim.t
        
        if σi == σR
            e = abs.(σ .- σi(model_i, load_sim))
        else
            e = abs.(σ .- σi(model, load_sim))
        end

        push!(errors, e)
        push!(timepoints, t)
    end 

    # Plotting and saving image
    plt = plot(timepoints, errors, label=["0.01" "0.03" "0.05" "0.07" "0.09"], 
                linewidth=2.0, legendfontsize = 10)
    xlabel!("time point")
    ylabel!("absolute error")
    title!("Changing step size. Deriv method = " * type[i])
    savefig("images/" * "deriv_eval_error_" * type[i] * ".svg")
    display(plt)

end



