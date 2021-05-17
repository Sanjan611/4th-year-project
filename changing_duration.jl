using RHEOS
using SpecialFunctions
using BenchmarkTools
using DataFrames
using CSV

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
σfs = [σA, σB, σR, σL1, σL12]
objs = [objA, objB, objR, objL1, objL12]

# DataFrame columns
sigma, other = [], []
duration = []
p0_c, p0_b, lo_c, lo_b, hi_c, hi_b = [], [], [], [], [], []
cs, bs, error, time_total = [], [], [], [] 

durations = range(5, 15, length = 10)

################## EXPERIMENT #################

for (k, dur) in enumerate(collect(durations))
    println("============= ",k," =============")

    time_sim = timeline(t_start = 0, t_end = dur, step = 0.05)
    load_sim = strainfunction(time_sim, t->f(t,model_params.β)); # f is the function from the paper

    for (i, σf) in enumerate(σfs) 
        println("============= ",k," ============= ",i, " ============= ")

        # Get the σ
        if σf!=σR
            σ = σf(model, load_sim)
        else
            σ = σf(model_i, load_sim)
        end

        data = RheoTimeData(σ = σ, ϵ = load_sim.ϵ, t = load_sim.t)

        for (j, obj) in enumerate(objs)
            println("============= ",k," ============= ",i, " ============= ", j, " ============= ")

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

            push!(sigma, label[i])
            push!(other, label[j])
            push!(duration, dur)
            push!(p0_c, p0.cᵦ)
            push!(p0_b, p0.β)
            push!(lo_c, lo.cᵦ)
            push!(lo_b, lo.β)
            push!(hi_c, hi.cᵦ)
            push!(hi_b, hi.β)
            push!(cs, minx[1])
            push!(bs, minx[2])
            push!(error, minf)
            push!(time_total, timetaken)
            
        end
    
    end
end

################## SAVE THE DATA TO A CSV FILE ##########################

df = DataFrame(sigma = sigma, other = other,
                duration = duration,
                p0_c = p0_c, p0_b = p0_b,
                lo_c = lo_c, lo_b = lo_b,
                hi_c = hi_c, hi_b = hi_b,
                cs = cs, bs = bs, error = error, time_total = time_total)

CSV.write("springpot_1_durations.csv", df)