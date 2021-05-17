using RHEOS
using SpecialFunctions
using BenchmarkTools
using DataFrames
using CSV
using Plots

include("optim.jl")
include("springpot_helper.jl")

############# OPTIMISATION INITIALISATION ##############
p0 = (cᵦ = 2.0, β = 0.5)
lo = (cᵦ = 0.0, β = 0.01)
hi = (cᵦ = 100.0, β = 0.99)

################### DEFINE SOME OTHER USEFUL VARIABLES #########################
label = ["A", "B", "R", "L1", "L12"] # the different methods used to evaluate
objs = [objA, objB, objR, objL1, objL12]

# arrays to store time and error
timetakens = [[], [], [], [], []]
errors = [[], [], [], [], []]
pred_βs = [[], [], [], [], []]


βs = range(0.01, 0.99, length = 30)

######### EXPERIMENT #################

for (k, β) in enumerate(collect(βs))
    println("============= ",k," =============")

    # Model and RheoTimeData 
    model_params = (cᵦ = 1.8, β = β)
    model = RheoModel(Springpot, model_params)
    model_i = RheoModel(Springpot_i, model_params)

    time_sim = timeline(t_start = 0, t_end = 8, step = 0.05)
    load_sim = strainfunction(time_sim, t->f(t,model_params.β)); # f is the function from the paper

    σ = σA(model, load_sim) # create simulated data using analytical expression
    data = RheoTimeData(σ = σ, ϵ = load_sim.ϵ, t = load_sim.t)

    for (j, obj) in enumerate(objs)
        println("============= ",k," ============= ",j, " ============= ")

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
        push!(errors[j], minf)
        push!(pred_βs[j], minx[2])
            
        
    
    end
end

###################### PLOTS ##############################

legend_labels = ["A" "B" "R" "L1" "L12"]

plt1 = plot(collect(βs), timetakens, ls = :auto, label = legend_labels, yaxis=:log, linewidth = 3, legend=:outertopright, legendfontsize = 10)
xlabel!("β")
ylabel!("Time taken for fitting")
# title!(" ")
savefig("images/A_timetaken_beta.svg")
display(plt1)

plt2 = plot(collect(βs), errors, ls = :auto, label = legend_labels, yaxis=:log, linewidth = 3, legend=:outertopright, legendfontsize = 10)
xlabel!("β")
ylabel!("Error from fit")
# title!(" ")
savefig("images/A_error_beta.svg")
display(plt2)

plt3 = plot(collect(βs), pred_βs, ls = :auto, label = legend_labels, linewidth = 3, legend=:outertopright, legendfontsize = 10)
xlabel!("β")
ylabel!("β after fitting")
# title!(" ")
savefig("images/A_predbeta_beta.svg")
display(plt3)

