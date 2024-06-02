using Revise
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons, SingleBumpAttractor
using Random, Distributions
using Q1
using BumpAttractorUtils
using Utils

np.tau = 3.0
N = np.N
n = sp.n
T = sp.T
delta_t = sp.delta_t

x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)
h_init = rand(Uniform(0,1), N) # initial potential values sampled from the uniform distribution

# plot parameters
nb_ticks_x = 5

# # plot I_ext for different times
# t_values = [350, 500, 650]
# p1 = plot()
# plot!(xlabel=L"x_i", ylabel=L"I_{ext}")
# for t in t_values
#     I_ext_t = Q1.I_ext.(x_i, t)
#     plot!(x_i, I_ext_t, label="t=$t", lw=2)
# end
# plot(p1)
# savefig("data/Q14_I_ext.pdf")

# simulate the network activity
spikes = SingleBumpAttractor.simulate_network(h_init, x_i, Q1.I_ext, 0.0, sp, np)

p = Utils.raster_plot(spikes, sp, np, title=LaTeXString("\$\\tau=$(np.tau)\$"))

map_time_to_idx(t) = Int(round(t/sp.delta_t))
plot!(map_time_to_idx.([300, 300, 400, 400]), [0, np.N, np.N, 0], lw=1, fill=(0, 0.25, :lightblue), color=:transparent, label=L"$I_\mathrm{ext}\neq 0$")
plot!(map_time_to_idx.([600, 600, 700, 700]), [0, np.N, np.N, 0], lw=1, fill=(0, 0.25, :lightblue), label=false, color=:transparent)

bin_size = 10 # ms, as advised in the instructions
Utils.plot_avg_bump_location(spikes, x_i, bin_size, sp, np)
# plot horizontal lines at 2*pi/3 and 4*pi/3 
idx = Utils.map_angle_to_idx(2*pi/3, np.N)
plot!([0, sp.n], [idx, idx], lw=1, ls=:dash, label=false, color=:black)
idx = Utils.map_angle_to_idx(4*pi/3, np.N)
plot!([0, sp.n], [idx, idx], lw=1, ls=:dash, label=false, color=:black)

# update y ticks
yticks!((
    Utils.map_angle_to_idx.([0, 2pi/3, pi, 4pi/3, 2pi], np.N),
    [L"0", L"\frac{2\pi}{3}", L"\pi", L"\frac{4\pi}{3}", L"2 \pi"])
)



display(p)
tau_str = replace(string(np.tau), "." => "_")
savefig("data/Q14_tau_$tau_str.pdf")
