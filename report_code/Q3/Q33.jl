using Revise
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Random, Distributions

using BumpAttractorUtils, SingleBumpAttractor, CoupledBumpAttractors
using Q3
using RandomTrajectory
using Utils

N = np.N
n = sp.n
T = sp.T
delta_t = sp.delta_t

# plot parameters
nb_ticks_x = 5
bin_size = 10 # ms, as advised in the instructions
bin_length = Int64(bin_size/delta_t)

sp.delay = 500
theta_coupled = deg2rad(10) # 10 degres

h_init_L = rand(Uniform(0,1), N) 
h_init_R = rand(Uniform(0,1), N)
h_init_H = rand(Uniform(0,1), N)
x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)

# define external input 
angles, pos = RandomTrajectory.random_trajectory(Q3.speed, T, delta_t, Q3.volatility)
plot(pos[:,1], pos[:,2], label=false)
Io = 1
t_to_idx(t) = Int(floor(t/delta_t)) + 1
theta(t) = angles[t_to_idx(t)]
I_ext(x, t) = Q3.I_ext_head(x, t, Io, theta)

# simulate the single bump attractor representing the head directions
S_i = SingleBumpAttractor.simulate_network(h_init_H, x_i, I_ext, 0.0, sp, np)
Utils.raster_plot(S_i, sp, np)

I_head_val = (S_i/delta_t) * cos.(x_i)
I_head(t) = I_head_val[t_to_idx(t)]

J_head = 1
I_ext_L(x, t) = -J_head/N * I_head(t)
I_ext_R(x, t) = J_head/N * I_head(t)

# coupled bump attractors
spikes_L, spikes_R = CoupledBumpAttractors.simulate_network(h_init_L, h_init_R, I_ext_L, I_ext_R, x_i, theta_coupled, sp, np)
spikes = (spikes_L + spikes_R)./2

heatmap(transpose(spikes), title="Network Activity", xlabel=L"t"*" (ms)", ylabel= "Neuron Location", c = reverse(cgrad(:grayC)), colorbar=false, right_margin = 3Plots.mm, left_margin = 2Plots.mm, yticks = (range(start = 0, stop = N , length =5), [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), xticks = (Int.(0:n/nb_ticks_x:n), Int.(0:T/nb_ticks_x:T)))
bump_location_L = locate_bump.(eachrow(spikes_L), Ref(x_i))
bump_location_bins_L = transpose(reshape(bump_location_L[1:n], bin_length, Int((n)/bin_length)))
avg_bump_location_L = locate_bump_avg.(Ref(ones(bin_length)), eachrow(bump_location_bins_L))
plot!(0:bin_length:n-1, avg_bump_location_L * (N/(2*pi)), label = "center of the bump")

bump_location_R = locate_bump.(eachrow(spikes_R), Ref(x_i))
bump_location_bins_R = transpose(reshape(bump_location_R[1:n], bin_length, Int((n)/bin_length)))
avg_bump_location_R = locate_bump_avg.(Ref(ones(bin_length)), eachrow(bump_location_bins_R))
plot!(0:bin_length:n-1, avg_bump_location_R * (N/(2*pi)), label = "center of the bump")