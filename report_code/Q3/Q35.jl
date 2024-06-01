using Revise
println(Revise.errors())
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Random, Distributions

using BumpAttractorUtils, SingleBumpAttractor, CoupledBumpAttractors
using Q3
using RandomTrajectory
using Utils

angles, pos = RandomTrajectory.random_trajectory(Q3.speed, sp.T, sp.delta_t, Q3.volatility)

x_i = collect(range(start = 0, stop = 2*pi, length = np.N + 1)[1:np.N]) # equally spaced neurons over the range [0, 2pi)
h_init = rand(Uniform(0,1), np.N) # initial potential values sampled from the uniform distribution

sp.delay = 500

t_to_idx(t) = Int(floor(t/sp.delta_t)) + 1
theta(t) = angles[t_to_idx(t)]
I_ext(x, t) = Q3.I_ext_head(x, t, Q3.Io, theta)

spikes = SingleBumpAttractor.simulate_network(h_init, x_i, angles, pos, I_ext, sp, np)

J_head = 2

I_head_x(t) = Q3.get_I_head_cos(spikes, t)
I_head_y(t) = Q3.get_I_head_sin(spikes, t)

spikes_L_x, spikes_R_x = CoupledBumpAttractors.integrate_position(J_head, I_head_x)
spikes_L_y, spikes_R_y = CoupledBumpAttractors.integrate_position(J_head, I_head_y)

spikes_x = (spikes_L_x + spikes_R_x)./2
spikes_y = (spikes_L_y + spikes_R_y)./2

bump_location_x = Utils.spikes_to_average_bump_location(spikes_x, x_i, Q3.bin_size, sp)
bump_location_y = Utils.spikes_to_average_bump_location(spikes_y, x_i, Q3.bin_size, sp)

plot(bump_location_x, bump_location_y)


