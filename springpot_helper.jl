using SpecialFunctions
using RHEOS
include("caputo.jl")

# STRAIN FUNCTION 
function f(t, k)
    # return t.^(k+β)
    return t.^(k)
end

# STRESS FUNCTION
function σA(model::RheoModel, load_sim::RheoTimeData, k)
    # return model.params.cᵦ*(gamma(k+model.params.β+1)/gamma(k+1))*(load_sim.t.^k)
    return model.params.cᵦ*(gamma(k+1)/gamma(k-model.params.β+1))*(load_sim.t.^(k-model.params.β))
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

function objA(data::RheoTimeData, model::RheoModel, params, g; verbose::Bool = false)
  
    time = copy(data.t)
    num_points = size(time,1)

    cᵦ = params[1]
    β = params[2]
    βₒ = model.params.β

    measured = data.σ
    k = 2
    ϵfdot = (gamma(k+1+βₒ)/gamma(k+1+βₒ-β))*(time.^(k+1+βₒ-β)) # OLD 
    cost = sum((measured .- cᵦ*ϵfdot).^2)/num_points

    
    if verbose
        println("cᵦ ", round(cᵦ, digits = 5), " ", "β ", round(β, digits = 5), " ", cost)
    end

    return cost
end

function objB(data::RheoTimeData, model::RheoModel, params, g, verbose::Bool = false)
    strain = copy(data.ϵ)
    time = copy(data.t)
    dt = time[2] - time[1]
    num_points = size(time,1)


    cᵦ = params[1]
    β = params[2]

    measured = data.σ
    ϵfdot = L1(strain, time, dt, β)
    
    cost = sum((measured .- cᵦ*ϵfdot).^2)/num_points

    if verbose
        println("cᵦ ", round(cᵦ, digits = 5), " ", "β ", round(β, digits = 5), " ", cost)
    end
    return cost
    return 0
end

function objR(model_i::RheoModel, load_sim::RheoTimeData)
    return 0
end

function objL1(data::RheoTimeData, model::RheoModel, params, g, verbose::Bool = false)
    strain = copy(data.ϵ)
    time = copy(data.t)
    dt = time[2] - time[1]
    num_points = size(time,1)


    cᵦ = params[1]
    β = params[2]

    measured = data.σ
    ϵfdot = L1(strain, time, dt, β)
    
    cost = sum((measured .- cᵦ*ϵfdot).^2)/num_points

    if verbose
        println("cᵦ ", round(cᵦ, digits = 5), " ", "β ", round(β, digits = 5), " ", cost)
    end
    return cost

end

function objL12(data::RheoTimeData, model::RheoModel, params, g, verbose::Bool = false)

    strain = copy(data.ϵ)
    time = copy(data.t)
    dt = time[2] - time[1]
    num_points = size(time,1)

    cᵦ = params[1]
    β = params[2]

    measured = data.σ
    ϵfdot = L12(strain, time, dt, β)
    cost = sum((measured .- cᵦ*ϵfdot).^2)/num_points

    if verbose
        println("cᵦ ", round(cᵦ, digits = 5), " ", "β ", round(β, digits = 5), " ", cost)
    end

    return cost
    
end




