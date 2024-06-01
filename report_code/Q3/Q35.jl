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

spikes = SingleBumpAttractor.simulate_network(h_init, x_i, I_ext, 0.0, sp, np)

J_head = 0.85

I_head_val_cos = (spikes/sp.delta_t) * cos.(x_i)
I_head_val_sin = (spikes/sp.delta_t) * sin.(x_i)
I_head_x(t) = I_head_val_cos[Q3.t_to_idx(t)]
I_head_y(t) = I_head_val_sin[Q3.t_to_idx(t)]

spikes_L_x, spikes_R_x = Q3.integrate_position(J_head, I_head_x)
spikes_L_y, spikes_R_y = Q3.integrate_position(J_head, I_head_y)

spikes_x = (spikes_L_x + spikes_R_x)./2
spikes_y = (spikes_L_y + spikes_R_y)./2

bin_size = 10
bump_location_x = Utils.spikes_to_average_bump_location(spikes_x, x_i, bin_size, sp) .- pi
bump_location_y = Utils.spikes_to_average_bump_location(spikes_y, x_i, bin_size, sp) .- pi

plot(bump_location_x, bump_location_y)
plot!(pos[:,1], pos[:,2])


