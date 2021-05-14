using SpecialFunctions
using RHEOS
include("caputo.jl")

# STRAIN FUNCTION 
function f(t, β)
    return t.^(1+β)
end

# STRESS FUNCTION
function σA(model::RheoModel, load_sim::RheoTimeData)
    return model.params.cᵦ*(gamma(2+model.params.β)/gamma(2))*load_sim.t
end

function σB(model::RheoModel, load_sim::RheoTimeData)
    data = modelpredict(load_sim, model, diff_method = "BD")
    return data.σ
end

function σR(model_i::RheoModel, load_sim::RheoTimeData)
    data = modelpredict(load_sim, model_i, diff_method = "BD")
    return data.σ
end

function σL1(model::RheoModel, load_sim::RheoTimeData)
    Δt = load_sim.t[2] - load_sim.t[1]
    σ = model.params.cᵦ*L1(load_sim.ϵ, load_sim.t, Δt, model.params.β)
    return σ
end

function σL12(model::RheoModel, load_sim::RheoTimeData)
    Δt = load_sim.t[2] - load_sim.t[1]
    σ = model.params.cᵦ*L12(load_sim.ϵ, load_sim.t, Δt, model.params.β)
    return σ
end

function objA(data::RheoTimeData, model::RheoModel, params, g)

    time = copy(data.t)
    num_points = size(time,1)

    cᵦ = params[1]
    β = params[2]
    βₒ = model.params.β

    ϵfdot = (gamma(2+βₒ)/gamma(2+βₒ-β))*(time.^(1+βₒ-β))
    cost = sum((measured .- cᵦ*ϵfdot).^2)/num_points

    return cost
end

function objB(data::RheoTimeData, model::RheoModel, params, g)

    data = modelpredict(data, model)
    return data.σ
end

function objR(model_i::RheoModel, load_sim::RheoTimeData)
    data = modelpredict(load_sim, model_i)
    return data.σ
end

function objL1(model::RheoModel, load_sim::RheoTimeData)
    strain = copy(load_sim.ϵ)
    time = copy(load_sim.t)
    dt = time[2] - time[1]

    ϵfdot = L1(strain, time, dt, model.params.β)
    return model.params.

end

function otherL12(model::RheoModel, load_sim::RheoTimeData)
    
end

function cost(data, params, g)
end


