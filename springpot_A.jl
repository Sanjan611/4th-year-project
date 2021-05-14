# A for analytical
# so σ(t) is obtained analytically

using RHEOS
using SpecialFunctions
include("optim.jl")

include("springpot_helper.jl")

# CREATE THE RHEOTIMEDATA

model_params = (cᵦ = 1.6, β = 0.5)
model = RheoModel(Springpot, model_params)
model_i = RheoModel(Springpot_i, model_params)

time_sim = timeline(t_start = 0, t_end = 8, step = 0.05)
load_sim = strainfunction(time_sim, t->f(t,model_params.β)); # f is the function from the paper

# t = data.t
# dt = t[2] - t[1]
# σ = data.σ 
# ϵ = data.ϵ

# σ = convert(Vector{Float64}, σ)
# ϵ = convert(Vector{Float64}, ϵ)
# t = convert(Vector{Float64}, t)

σfs = [σA, σB, σR, σL1, σL12]
others = [otherA, otherB, otherR, otherL1, otherL12]


for (i, σf) in enumerate(σfs) 

    # Get the σ
    if σf!=σR
        σ = σf(model, copy(load_sim))
        print("True")
    else
        σ = σf(model_i, copy(load_sim))
    end

    data = RheoTimeData(σ = σ, ϵ = load_sim.ϵ, t = load_sim.t)

    for (j, other) in enumerate(others)
        obj = 
        myleastsquares(params_init = [2.0, 0.7], low_bounds = [0.0, 0.0], hi_bounds = [100, 1.0], data = data, obj = )
    end
 


  
end

