using Revise
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

Io = 0
t_to_idx(t) = Int(floor(t/delta_t)) + 1
theta(t) = angles[t_to_idx(t)]
I_ext(x, t) = Q3.I_ext_head(x, t, Io, theta)

S_i = SingleBumpAttractor.simulate_network(x_i, h_init, I_ext, 0.0, sp, np)

plot(angles, label=false)
angle_segments = Utils.split_into_segments(angles)
plot(angle_segments, label=false)
count = 1
p = plot()
for segment in angle_segments
    new_n = count + length(segment)-1
    plot!(n:new_n, segment, label=false, plot=p, color=:black, lw=1)
    count = new_n
end
display(p)
# savefig("data/Q32_angles.png")

# find the location of the bump at every timestep
bump_location = locate_bump.(eachrow(S_i), Ref(x_i))
# average location of the bump over small bins of time (to have a clearer plot)
bin_size = 10 # ms, as advised in the instructions
bin_length = Int64(bin_size/delta_t)
bump_location_bins = transpose(reshape(bump_location[1:n], bin_length, Int((n)/bin_length)))
avg_bump_location = locate_bump_avg.(Ref(ones(bin_length)), eachrow(bump_location_bins)) # need to use a circular mean method again
Utils.raster_plot(S_i, sp, np)
plot!(0:bin_length:n-1, avg_bump_location * (N/(2*pi)), label = "center of the bump")

