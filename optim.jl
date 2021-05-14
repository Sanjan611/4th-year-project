using NLopt
using BenchmarkTools

function myleastsquares(;params_init::Array{Float64, 1}, low_bounds::Array{Float64, 1}, hi_bounds::Array{Float64, 1},
    data::RheoTimeData, model::RheoModel, obj::Function, _rel_tol = 1e-4)

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

    min_objective!(opt, (params_init, g) -> obj(data, model, params_init, g))

    (minf, minx, ret) = NLopt.optimize(opt, params_init)

    numevals = opt.numevals # the number of function evaluations
    println("got $minf at $minx after $numevals iterations (returned $ret)")

    return (minf, minx, ret)

end
