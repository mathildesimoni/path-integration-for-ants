using Revise
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using CoupledBumpAttractors, Neurons, BumpAttractorUtils
using Random, Distributions
using Q2
using Utils

N = np.N
n = sp.n
T = sp.T
x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)
I_ext(x::Real, t::Real) = 0.0 # no external input

# plot parameters
nb_ticks_x = 5
bin_size = 10 # ms, as advised in the instructions
bin_length = Int64(bin_size/sp.delta_t)

# coupled bump attractor parameters
theta = deg2rad(10) # 10 degres
h_init_L = rand(Uniform(0,1), N) 
h_init_R = rand(Uniform(0,1), N)

# simulate the network activity
spikes_L, spikes_R = CoupledBumpAttractors.simulate_network(h_init_L, h_init_R, I_ext, I_ext, x_i, theta, sp, np)
avg_bump_location_L = Utils.spikes_to_average_bump_location(spikes_L, x_i, bin_size, sp)
avg_bump_location_R = Utils.spikes_to_average_bump_location(spikes_R, x_i, bin_size, sp)

# plot
spikes = (spikes_L + spikes_R)./2 # for the plots
p = Utils.raster_plot(spikes, sp, np)
Utils.plot_segments(collect(0:bin_length:n-1), avg_bump_location_L * (N/(2*pi)), color = :blue, label = "center of the left bump")
Utils.plot_segments(collect(0:bin_length:n-1), avg_bump_location_R * (N/(2*pi)), color = :red, label = "center of the right bump")
display(p)
# savefig("data/Q22.pdf")