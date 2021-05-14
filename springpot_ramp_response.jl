using RHEOS
using Plots

Plots.pyplot()

model_params = (cᵦ = 1.5, β=0.5)
model = RheoModel(Springpot, model_params)
time_sim = timeline(t_start = 0, t_end = 8, step = 0.1);
load_sim = strainfunction(time_sim, ramp(offset = 2.0, gradient = 3.0));
# load_sim = strainfunction(time_sim, hstep(offset = 2.0, amp = 1.0));

data = modelpredict(load_sim, model, diff_method = "BD")

plot(data.t, data.σ)
# plot!(data.t, data.ϵ)
# legend(["data", "data"])