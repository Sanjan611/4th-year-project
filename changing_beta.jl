using RHEOS
using SpecialFunctions
using BenchmarkTools
using Plots

include("optim.jl")
include("springpot_helper.jl")

############# OPTIMISATION INITIALISATION ##############
p0 = (cᵦ = 2.0, β = 0.5)
lo = (cᵦ = 0.0, β = 0.01)
hi = (cᵦ = 10.0, β = 0.99)

################### DEFINE SOME OTHER USEFUL VARIABLES #########################
# label = ["A", "B", "R", "L1", "L12"] # the different methods used to evaluate
# objs = [objA, objB, objR, objL1, objL12]
label = ["B", "R", "L1", "L12"] # the different methods used to evaluate
objs = [objB, objR, objL1, objL12]

# arrays to store time and error
timetakens = [[], [], [], []]
errors = [[], [], [], []]
pred_βs = [[], [], [], []]


βs = range(0.01, 0.99, length = 20)

######### EXPERIMENT #################
λ = 2.5
for (i, β) in enumerate(collect(βs))
    println("============= ",i," =============")

    # Model and RheoTimeData 
    model_params = (cᵦ = 1.8, β = β)
    model = RheoModel(Springpot, model_params)
    model_i = RheoModel(Springpot_i, model_params)

    time_sim = timeline(t_start = 0, t_end = 8, step = 0.05)
    load_sim = strainfunction(time_sim, t->f(t, λ)); # f is the function from the paper
    σ = σA(model, load_sim, λ) # create simulated data using analytical expression
    data = RheoTimeData(σ = σ, ϵ = load_sim.ϵ, t = load_sim.t)

    for (j, obj) in enumerate(objs)
        println("============= ",i," ============= ",j, " ============= ")

        if obj == objB
            fitted_model, timetaken, ret, minx, minf, numevals = 
                                    modelfit(data, Springpot, strain_imposed,
                                    p0 = p0,
                                    lo = lo, 
                                    hi = hi,
                                    return_stats = true,
                                    )
        
        elseif obj == objR
            fitted_model, timetaken, ret, minx, minf, numevals = 
                                    modelfit(data, Springpot_i, strain_imposed,
                                    p0 = p0,
                                    lo = lo, 
                                    hi = hi,
                                    return_stats = true
                                    )
        else
            (minf, minx, ret), timetaken, bytes, gctime, memalloc = 
                            @timed myleastsquares(params_init = [p0.cᵦ, p0.β], 
                            low_bounds = [lo.cᵦ, lo.β], hi_bounds = [hi.cᵦ, hi.β], 
                            data = data, model = model, obj = obj, verbose = false)
            println("Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")
        end

        push!(timetakens[j], timetaken)
        push!(errors[j], minf/σ[end])
        push!(pred_βs[j], minx[2])
            
        
    
    end
end

###################### PLOTS ##############################

legend_labels = ["B" "R" "L1" "L12"]

plt1 = plot(collect(βs), timetakens, ls = :auto, label = legend_labels, yaxis=:log, linewidth = 3, legend=:outertopright, legendfontsize = 10)
xlabel!("β")
ylabel!("(Log) Time taken (s)")
# title!(" ")
title!("λ = " * string(λ))
savefig("images/A_timetaken_beta_lambda_"*string(λ)*".svg")
display(plt1)

plt2 = plot(collect(βs), errors, ls = :auto, label = legend_labels, yaxis=:log, linewidth = 3, legend=:outertopright, legendfontsize = 10)
xlabel!("β")
ylabel!("(Log) Relative error")
# title!(" ")
title!("λ = " * string(λ))
savefig("images/A_error_beta_lambda_"*string(λ)*".svg")
display(plt2)

plt3 = plot(collect(βs), pred_βs, ls = :auto, label = legend_labels, linewidth = 3, legend=:outertopright, legendfontsize = 10)
xlabel!("β")
ylabel!("β after fitting")
# title!(" ")
title!("λ = " * string(λ))
savefig("images/A_predbeta_beta_lambda_"*string(λ)*".svg")
display(plt3)

