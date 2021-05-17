using RHEOS
using SpecialFunctions
using BenchmarkTools
using DataFrames
using CSV
using Plots

include("optim.jl")
include("springpot_helper.jl")

############# SET UP THE MODEL AND OPTIMISATION INITIALISATION ##############

# Model and RheoTimeData 
model_params = (cᵦ = 1.8, β = 0.3)
model = RheoModel(Springpot, model_params)
model_i = RheoModel(Springpot_i, model_params)

p0 = (cᵦ = 2.0, β = 0.7)
lo = (cᵦ = 0.0, β = 0.01)
hi = (cᵦ = 100.0, β = 0.99)

################### DEFINE SOME OTHER USEFUL VARIABLES #########################

label = ["A", "B", "R", "L1", "L12"] # the different methods used to evaluate
objs = [objA, objB, objR, objL1, objL12]

# arrays to store time and error
timetakens = [[], [], [], [], []]
errors = [[], [], [], [], []]

durations = range(5, 20, length = 20)


for (k, dur) in enumerate(collect(durations))
    println("============= ",k," =============")

    time_sim = timeline(t_start = 0, t_end = dur, step = 0.05)
    load_sim = strainfunction(time_sim, t->f(t,model_params.β)); # f is the function from the paper

    σ = σA(model, load_sim)
    data = RheoTimeData(σ = σ, ϵ = load_sim.ϵ, t = load_sim.t)

    for (j, obj) in enumerate(objs)
        println("============= ",k," ============= ",j, " ============= ")

        if obj == objB
            fitted_model, timetaken, ret, minx, minf, numevals = 
                                    modelfit(data, Springpot, strain_imposed,
                                    p0 = p0,
                                    lo = lo, 
                                    hi = hi,
                                    return_stats = true
                                    )
            println("Numevals: $numevals")
        
        elseif obj == objR
            fitted_model, timetaken, ret, minx, minf, numevals = 
                                    modelfit(data, Springpot_i, strain_imposed,
                                    p0 = p0,
                                    lo = lo, 
                                    hi = hi,
                                    return_stats = true
                                    )
            println("Numevals: $numevals")

        else
            (minf, minx, ret, numevals), timetaken, bytes, gctime, memalloc = 
                            @timed myleastsquares(params_init = [p0.cᵦ, p0.β], 
                            low_bounds = [lo.cᵦ, lo.β], hi_bounds = [hi.cᵦ, hi.β], 
                            data = data, model = model, obj = obj, verbose = false)
            println("Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf, Numevals: $numevals")
        end

        push!(timetakens[j], timetaken)
        push!(errors[j], minf)
            
        
    
    end
end

##################### PLOTS ########################

legend_labels = ["A" "B" "R" "L1" "L12"]

plt1 = plot(collect(durations), timetakens, ls = :auto, label = legend_labels, yaxis=:log, linewidth = 3, legend=:outertopright, legendfontsize = 10)
xlabel!("Experiment duration")
ylabel!("Time taken for fitting")
# title!("Time taken vs experiment duration")
savefig("images/A_timetaken_duration.svg")
display(plt1)

plt2 = plot(collect(durations), errors, ls = :auto, label = legend_labels, yaxis=:log, linewidth = 3, legend=:outertopright, legendfontsize = 10)
xlabel!("Experiment duration")
ylabel!("Error from fit")
# title!("Obj error vs experiment duration")
savefig("images/A_error_duration.svg")
display(plt2)

