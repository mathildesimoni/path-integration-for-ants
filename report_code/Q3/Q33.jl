using Revise
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Random, Distributions

using BumpAttractorUtils, SingleBumpAttractor, CoupledBumpAttractors
using Q3
using RandomTrajectory
using Utils

println(Revise.errors())

N = np.N
n = sp.n
T = sp.T
delta_t = sp.delta_t
x_i = Q3.x_i

# plot parameters
nb_ticks_x = 5
bin_size = 10 # ms, as advised in the instructions
bin_length = Int64(bin_size/delta_t)

sp.delay = 500

h_init_H = rand(Uniform(0,1), N)

# define external input 
angles, pos = RandomTrajectory.random_trajectory(Q3.speed, T, delta_t, Q3.volatility)
plot(pos[:,1], pos[:,2], label=false)
Io = 1
t_to_idx(t) = Int(floor(t/delta_t)) + 1
theta(t) = angles[t_to_idx(t)]
I_ext(x, t) = Q3.I_ext_head(x, t, Io, theta)

# simulate the single bump attractor representing the head directions
S_i = SingleBumpAttractor.simulate_network(h_init_H, x_i, I_ext, 0.0, sp, np)

I_head_val = (S_i/delta_t) * cos.(x_i)
I_head(t) = I_head_val[Q3.t_to_idx(t)]

J_head = 3

spikes_L, spikes_R = Q3.integrate_position(J_head, I_head)

spikes = (spikes_L + spikes_R)./2
p = Utils.raster_plot(spikes, sp, np)
Utils.plot_avg_bump_location(spikes_L, x_i, bin_size, sp, np, color=:red, label=LaTeXString(L"\theta_{\mathrm{bump}}^{\mathrm{L}}"))
Utils.plot_avg_bump_location(spikes_R, x_i, bin_size, sp, np, color=:blue, label=LaTeXString(L"\theta_{\mathrm{bump}}^{\mathrm{R}}"))
display(p)

plot(pos[:,1])

J_head_range = collect(range(start = -4, stop = 4, length = 51))
bump_locations = zeros(length(J_head_range))
sp.rate_neurons = true
for (i, J) in enumerate(J_head_range)
    spikes_L, spikes_R = Q3.integrate_position(J, I_head)
    spikes = (spikes_L + spikes_R)./2
    avg_bump_location_L = Utils.spikes_to_average_bump_location(spikes_L, x_i, bin_size, sp)
    avg_bump_location_R = Utils.spikes_to_average_bump_location(spikes_R, x_i, bin_size, sp)
    last_bump_location = 0.5*(avg_bump_location_L[end] + avg_bump_location_R[end])
    first_bump_location = 0.5*(avg_bump_location_L[1] + avg_bump_location_R[1])
    bump_locations[i] = last_bump_location - first_bump_location
end 

p = plot()
Utils.plot_segments(J_head_range, bump_locations)
display(p)

# plot(pos[:,1])