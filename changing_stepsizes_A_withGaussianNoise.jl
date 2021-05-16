using RHEOS
using SpecialFunctions
using BenchmarkTools
using DataFrames
using CSV
using Plots
using Noise

include("optim.jl")
include("springpot_helper.jl")

p0 = (cᵦ = 2.0, β = 0.7)
lo = (cᵦ = 0.0, β = 0.01)
hi = (cᵦ = 100.0, β = 0.99)

label = ["A", "B", "R", "L1", "L12"] # the different methods used to evaluate
objs = [objA, objB, objR, objL1, objL12]

# DataFrame columns
sigma, other = [], []
stepsize = []
p0_c, p0_b, lo_c, lo_b, hi_c, hi_b = [], [], [], [], [], []
cs, bs, error, time_total = [], [], [], [] 

stepsizes = range(0.05, 0.2, length = 30)

# Model and RheoTimeData 
model_params = (cᵦ = 1.8, β = 0.3)
model = RheoModel(Springpot, model_params)
model_i = RheoModel(Springpot_i, model_params)

# arrays to store time and error
timetakens = [[], [], [], [], []]
errors = [[], [], [], [], []]

for (k, ss) in enumerate(collect(stepsizes))
    println("============= ",k," =============")

    time_sim = timeline(t_start = 0, t_end = 10, step = ss)
    load_sim = strainfunction(time_sim, t->f(t,model_params.β)); # f is the function from the paper

    σ = σA(model, load_sim)
    data = RheoTimeData(σ = add_gauss(σ, 0.5), ϵ = load_sim.ϵ, t = load_sim.t)

    for (j, obj) in enumerate(objs)
        println("============= ",k," ============= ",j, " ============= ")

        if obj == objB
            fitted_model, timetaken, ret, minx, minf = 
                                    modelfit(data, Springpot, strain_imposed,
                                    p0 = p0,
                                    lo = lo, 
                                    hi = hi,
                                    return_stats = true,
                                    )
        
        elseif obj == objR
            fitted_model, timetaken, ret, minx, minf = 
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

        # push!(sigma, label[i])
        # push!(other, label[j])
        # push!(stepsize, ss)
        # push!(p0_c, p0.cᵦ)
        # push!(p0_b, p0.β)
        # push!(lo_c, lo.cᵦ)
        # push!(lo_b, lo.β)
        # push!(hi_c, hi.cᵦ)
        # push!(hi_b, hi.β)
        # push!(cs, minx[1])
        # push!(bs, minx[2])
        # push!(error, minf)
        # push!(time_total, timetaken)

        push!(timetakens[j], timetaken)
        push!(errors[j], minf)
            
        
    
    end
end

# df = DataFrame(sigma = sigma, other = other,
#                 stepsize = stepsize,
#                 p0_c = p0_c, p0_b = p0_b,
#                 lo_c = lo_c, lo_b = lo_b,
#                 hi_c = hi_c, hi_b = hi_b,
#                 cs = cs, bs = bs, error = error, time_total = time_total)

# CSV.write("springpot_1_stepsizes.csv", df)

legend_labels = ["A" "B" "R" "L1" "L12"]

plt1 = plot(collect(stepsizes), timetakens, ls = :auto, label = legend_labels, yaxis=:log, linewidth = 2.5, legend=:outertopright)
xlabel!("stepsize")
ylabel!("time taken")
title!("Time taken vs step size")
display(plt1)

plt2 = plot(collect(stepsizes), errors, ls = :auto, label = legend_labels, yaxis=:log, linewidth = 2.5, legend=:outertopright)
xlabel!("stepsize")
ylabel!("objective function error")
title!("Objective error vs step size")
display(plt2)

