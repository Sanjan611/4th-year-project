using Plots
# using Pkg
# Pkg.build("RHEOS")

using RHEOS
include("caputo.jl")
using SpecialFunctions

Plots.pyplot()

cᵦs = 1:0.04:3
βs = 0.4:0.01:0.8

# STRAIN FUNCTION 
function f(t)
    return t.^1.5
end

# STRESS FUNCTION
function fσ(t, cᵦ, β)
    return cᵦ*(gamma(2+β)/gamma(2))*t
end

# CREATE THE RHEOTIMEDATA
model_params = (cᵦ = 1.6, β = 0.5)
model = RheoModel(Springpot, model_params)
time_sim = timeline(t_start = 0, t_end = 8, step = 0.1)
load_sim = strainfunction(time_sim, f)

# data = modelpredict(load_sim, model, diff_method = "BD")

data = RheoTimeData(σ = fσ(load_sim.t, model_params.cᵦ, model_params.β), ϵ = load_sim.ϵ, t = load_sim.t)



σ = data.σ
ϵ = data.ϵ
t = data.t

count = 0

function Ga(m::RheoModel)
    if m.expressions.Ga_safe
        m._Ga
    else
        ta -> [m._G(t) for t in ta]
    end
end

relaxmodulus = Ga(model)

function springpot_cost_convint(cᵦ, β)

    params = (cᵦ= cᵦ, β = β)
    time_series = copy(t)
    println(params)

    modulus = relaxmod(Springpot, params)
    mod = modulus

    deriv = RHEOS.derivBD
    dcontrolled = deriv(ϵ, t)
    
    measured = copy(data.σ)

    dt = time_series[2] - time_series[1]
    println("dt = ", dt)

    # Singularity related check
    modsing = (t->model._G(t,params))
    sing = RHEOS.singularitytest(modsing)
    println("singularity: ", sing)
    if sing 
        time_series[1] = 0.0 + (time_series[2] - time_series[1])/10
    end

    convolved = RHEOS.boltzconvolve(mod, time_series, dt, dcontrolled)
    # println("convolved sum = ", sum(convolved))
    # println(convolved)
    num_points = length(convolved)
    println("num points = ", num_points)

    cost = sum((measured - convolved).^2)/num_points
    println("cost = ", cost)
    return cost
end

function springpot_cost_DE(cᵦ, β)

    global count +=  1 # just for printing purposes

    measured = copy(data.σ)
    time = copy(data.t)
    strain = copy(data.ϵ)


    cost = 0
    
    num_points = size(measured,1)
    dt = time[2] - time[1]

    # HOW TO CALCULATE THE FRACTIONAL DERIVATIVE BIT
    ϵfdot = L1(strain, time, dt, β)
    # βₒ = 0.5
    # ϵfdot = (gamma(2+βₒ)/gamma(2+βₒ-β))*(time.^(1+βₒ-β))

    # Cost is the mean squared error 
    cost = sum((measured .- cᵦ*ϵfdot).^2)/num_points

    println("β ", β, " ", "cᵦ ", cᵦ, " ", count, " ", cost)
    
    return cost
end

# CREATING THE CONTOUR MAP
X = repeat(reshape(cᵦs, 1, :), length(βs),1)
Y = repeat(βs, 1, length(cᵦs))

println("Now creating Z")
# Z = map(springpot_cost_convint, X, Y)
Z = map(springpot_cost_DE, X, Y)

p1 = contour(cᵦs, βs, Z, levels = collect(0:0.02:0.8), fill=true,clims=(0.,0.8), xlabel = "cᵦ", ylabel = "β")
scatter!(values(model_params))
# p2 = contour(cᵦ, β, Z, xlabel = "cᵦ", ylabel = "β")

plot(p1)

