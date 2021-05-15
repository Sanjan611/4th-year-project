# A for analytical
# so σ(t) is obtained analytically

using RHEOS
using SpecialFunctions
using BenchmarkTools

include("optim.jl")

include("springpot_helper.jl")

# CREATE THE RHEOTIMEDATA

model_params = (cᵦ = 1.6, β = 0.5)
model = RheoModel(Springpot, model_params)
model_i = RheoModel(Springpot_i, model_params)

time_sim = timeline(t_start = 0, t_end = 8, step = 0.05)
load_sim = strainfunction(time_sim, t->f(t,model_params.β)); # f is the function from the paper

σfs = [σA, σB, σR, σL1, σL12]
objs = [objA, objB, objR, objL1, objL12]


for (i, σf) in enumerate(σfs) 
    println("==================== ",i," ===================")

    # Get the σ
    if σf!=σR
        σ = σf(model, load_sim)
    else
        σ = σf(model_i, load_sim)
    end

    data = RheoTimeData(σ = σ, ϵ = load_sim.ϵ, t = load_sim.t)

    for (j, obj) in enumerate(objs)
        println("==================== ",i," =================== ", j, " =================== ")

        if obj == objB
            fitted_model = modelfit(data, Springpot, strain_imposed,
                                    p0 = (cᵦ = 2.0, β = 0.7),
                                    lo = (cᵦ = 0.0, β = 0.00), 
                                    hi = (cᵦ = 100.0, β = 1.0)
                                    )
        
        elseif obj == objR
            fitted_model = modelfit(data, Springpot_i, strain_imposed,
                                    p0 = (cᵦ = 2.0, β = 0.7),
                                    lo = (cᵦ = 0.0, β = 0.00), 
                                    hi = (cᵦ = 100.0, β = 1.0)
                                    )
        else
            (minf, minx, ret), timetaken, bytes, gctime, memalloc = 
                            @timed myleastsquares(params_init = [2.0, 0.7], 
                            low_bounds = [0.0, 0.0], hi_bounds = [100.0, 1.0], 
                            data = data, model = model, obj = obj)
            println("Time: $timetaken s, Why: $ret, Parameters: $minx, Error: $minf")
        end
        
        
        
    end
 


  
end

