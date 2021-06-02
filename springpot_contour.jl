using Plots
using RHEOS
using Noise
include("caputo.jl")
include("springpot_helper.jl")
using SpecialFunctions

cᵦs = 1:0.05:2
βs = 0.4:0.05:0.8

# CREATE THE RHEOTIMEDATA
model_params = (cᵦ = 1.6, β = 0.5)
model = RheoModel(Springpot, model_params)
model_i = RheoModel(Springpot_i, model_params)
time_sim = timeline(t_start = 0, t_end = 8, step = 0.1)
k = 1
load_sim = strainfunction(time_sim, t->f(t,model_params.β, k=k));

σ = σA(model, load_sim, k=k)
# data = RheoTimeData(σ = add_gauss(σ, 0.5), ϵ = load_sim.ϵ, t = load_sim.t)
data = RheoTimeData(σ = σ, ϵ = load_sim.ϵ, t = load_sim.t)

count = 0
function cost(cᵦ, β)
    println(count)
    global count = count + 1;

    model = RheoModel(Springpot, (cᵦ = cᵦ, β = β))
    err = (σ .- σA(model, load_sim))
    num_points = length(σ)
    mse = sum(err.^2) / num_points
    return mse
end


# CREATING THE CONTOUR MAP
X = repeat(reshape(cᵦs, 1, :), length(βs),1)
Y = repeat(βs, 1, length(cᵦs))

println("Now creating Z")
Z = map(cost, X, Y)

p1 = contour(cᵦs, βs, Z, levels = collect(0:0.02:0.8), fill=true,clims=(0.,0.8), xlabel = "cᵦ", ylabel = "β")
scatter!(values(model_params))
# p2 = contour(cᵦ, β, Z, xlabel = "cᵦ", ylabel = "β")

plot(p1)

