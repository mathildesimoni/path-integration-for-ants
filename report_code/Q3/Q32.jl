using Revise
println(Revise.errors())
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Distributions

using BumpAttractorUtils, SingleBumpAttractor
using Q3
using RandomTrajectory
using Utils

N = np.N
n = sp.n
T = sp.T
delta_t = sp.delta_t
bin_size = 10 # ms, as advised in the instructions
bin_length = Int64(bin_size/sp.delta_t)

angles, pos = RandomTrajectory.random_trajectory(Q3.speed, T, delta_t, Q3.volatility)

x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)
h_init = rand(Uniform(0,1), N) # initial potential values sampled from the uniform distribution

Io = 1
t_to_idx(t) = Int(floor(t/delta_t)) + 1
theta(t) = angles[t_to_idx(t)]
I_ext(x, t) = Q3.I_ext_head(x, t, Io, theta)

S_i = SingleBumpAttractor.simulate_network(h_init, x_i, I_ext, 0.0, sp, np)

p = Utils.raster_plot(S_i, sp, np)
bin_size = 5 # ms, as advised in the instructions
Utils.plot_angle_location(angles, N, p, color = :lightblue, label=LaTeXString(L"\theta^{\mathrm{H}}"))
Utils.plot_avg_bump_location(S_i, x_i, bin_size, sp, np, color = :red, label=LaTeXString(L"\theta_{\mathrm{bump}}^{\mathrm{H}}"))
display(p)

