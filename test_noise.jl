using RHEOS
using PyPlot

PyPlot.clf()

# Step 1: Create a model and time data
# model = RheoModel(KelvinVoigt, k=1, η=1)
model = RheoModel(Maxwell, k=1, η=1)
data_t = timeline(t_start=0, t_end=20.0, step = 0.02)

# Step 2: Create strain data (noise free)
# data = strainfunction(data_t,hstep(offset=5.0,amp=10));
data = strainfunction(data_t,ramp(offset=2.0,gradient=0.8));

# Step 3: Get stress response
data_p = modelpredict(data, model)

# Step 4: Create same strain data with added noise
noise = strainfunction(data_t, t -> rand());
data_noisy = data + noise
data_smooth = smooth(data_noisy, 1);

# Step 5: Get stress response (with noise)
data_p_noisy = modelpredict(data_noisy, model);
data_p_smooth = modelpredict(data_smooth, model);

# Step 6: Quantify the noise (output/input)
plot(data_p.t, data_p.σ)
plot(data_p_smooth.t, data_p_smooth.σ)

# plot(data.t, data.ϵ)
# plot(data_smooth.t, data_smooth.ϵ)

legend(["True stress response", "Noisy stress response"])

xlabel("Time")
ylabel("Stress")

savefig("images/test_noise.eps")

display(gcf())

