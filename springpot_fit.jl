using RHEOS
using NLopt
using BenchmarkTools

include("caputo.jl")

function springpot_cost_DE(data, params, g)


    cᵦ = params[1]
    β = params[2]

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

    # Cost is the mean squared error # Diff eq
    cost = sum((measured .- cᵦ*ϵfdot).^2)/num_points

    println("cᵦ ", round(cᵦ, digits = 5), " ", "β ", round(β, digits = 5), " ", count, " ", cost)
    
    return cost
end

function myleastsquares(;params_init::Array{Float64, 1}, low_bounds::Array{Float64, 1}, hi_bounds::Array{Float64, 1},
    data::RheoTimeData, _rel_tol = 1e-4)

    # initialise NLOpt.Opt object with :LN_SBPLX Subplex algorithm
    opt = Opt(:LN_SBPLX, length(params_init))

    # set lower bounds and upper bounds unless they take null value
    if !isnothing(low_bounds)
        low_bounds = convert(Vector{Float64},low_bounds)
        lower_bounds!(opt, low_bounds)
    end

    if !isnothing(hi_bounds)
        hi_bounds = convert(Vector{Float64}, hi_bounds)
        upper_bounds!(opt, hi_bounds)
    end

    # set relative tolerance
    xtol_rel!(opt, _rel_tol)

    println(opt.xtol_rel)

    # Convert to float64 to avoid conversion by NLOpt
    params_init = convert(Vector{Float64},params_init)

    min_objective!(opt, (params_init, g) -> springpot_cost_DE(data, params_init, g))

    # println(opt.min_objective) # no readable property like this

    (minf, minx, ret) = NLopt.optimize(opt, params_init)

    numevals = opt.numevals # the number of function evaluations
    println("got $minf at $minx after $numevals iterations (returned $ret)")

    return (minf, minx, ret)

end

# STRAIN FUNCTION 
function f(t)
    return t.^1.5
end

count = 0

# CREATE THE RHEOTIMEDATA
model_params = (cᵦ = 1.6, β = 0.3)
model = RheoModel(Springpot, model_params)
time_sim = timeline(t_start = 0, t_end = 8, step = 0.05)
# load_sim = strainfunction(time_sim, f)
# load_sim = strainfunction(time_sim, ramp(offset = 0.0, gradient = 1.0));
load_sim = strainfunction(time_sim, hstep(offset = 1.0, amp = 1.0));

data = modelpredict(load_sim, model, diff_method = "BD")

t = data.t
dt = t[2] - t[1]
σ = data.σ 
ϵ = data.ϵ

σ = convert(Vector{Float64}, σ)
ϵ = convert(Vector{Float64}, ϵ)
t = convert(Vector{Float64}, t)




myleastsquares(params_init = [2.0, 0.7], low_bounds = [0.0, 0.0], hi_bounds = [100, 1.0], data = data)