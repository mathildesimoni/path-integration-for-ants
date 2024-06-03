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
volatility = Q3.volatility 
volatility = 0.08
x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)

# plot parameters
nb_ticks_x = 5
bin_size = 10 # ms, as advised in the instructions
bin_length = Int64(bin_size/delta_t)

sp.delay = 500

h_init_H = rand(Uniform(0,1), N)

num_trajectories = 100
max_I_heads = zeros(num_trajectories)
min_I_heads = zeros(num_trajectories)
max_val = 0
min_val = 0
traces = zeros(num_trajectories, n+1)
for i in 1:num_trajectories
    # define external input 
    angles, pos = RandomTrajectory.random_trajectory(Q3.speed, T, delta_t, volatility)
    # RandomTrajectory.plot_trajectory(pos[:,1], pos[:,2])
    # p = plot!()
    # display(p)
    
    Io = 1
    t_to_idx(t) = Int(floor(t/delta_t)) + 1
    theta(t) = angles[t_to_idx(t)]
    I_ext(x, t) = Q3.I_ext_head(x, t, Io, theta)
    
    # simulate the single bump attractor representing the head directions
    S_i = SingleBumpAttractor.simulate_network(h_init_H, x_i, I_ext, 0.0, sp, np)
    
    I_head_val = (S_i./delta_t) * cos.(x_i)
    
    trace = I_head_val./np.N
    
    max_I_heads[i] = maximum(trace)
    min_I_heads[i] = minimum(trace)
    traces[i, :] = trace
    max_val = max(max_val, max_I_heads[i])
    min_val = min(min_val, min_I_heads[i])
    println("Trajectory $i,\t max: $(max_I_heads[i]),\t min: $(min_I_heads[i])")
end 

plot()
xlabel!("Trajectory")
ylabel!(L"I^R_{\mathrm{head}}")
# hline!([max_val], color=:red, label="Max", linestyle=:dash)
# hline!([min_val], color=:blue, label="Min", linestyle=:dash)
scatter!(max_I_heads, label="maximum", color=:red)
scatter!(min_I_heads, label="minimum", color=:blue)
savefig("data/Q33_I_head.pdf")

maximum(max_I_heads), minimum(min_I_heads)

J_head = 1

spikes_L, spikes_R = Q3.integrate_position(J_head, I_head)

spikes = (spikes_L + spikes_R)./2
p = Utils.raster_plot(spikes, sp, np)
Utils.plot_avg_bump_location(spikes_L, x_i, bin_size, sp, np, color=:red, label=LaTeXString(L"\theta_{\mathrm{bump}}^{\mathrm{L}}"))
Utils.plot_avg_bump_location(spikes_R, x_i, bin_size, sp, np, color=:blue, label=LaTeXString(L"\theta_{\mathrm{bump}}^{\mathrm{R}}"))
display(p)

# J_head_range = collect(range(start = -4, stop = 4, length = 51))
# bump_locations = zeros(length(J_head_range))
# sp.rate_neurons = true
# for (i, J) in enumerate(J_head_range)
#     spikes_L, spikes_R = Q3.integrate_position(J, I_head)
#     spikes = (spikes_L + spikes_R)./2
#     avg_bump_location_L = Utils.spikes_to_average_bump_location(spikes_L, x_i, bin_size, sp)
#     avg_bump_location_R = Utils.spikes_to_average_bump_location(spikes_R, x_i, bin_size, sp)
#     last_bump_location = 0.5*(avg_bump_location_L[end] + avg_bump_location_R[end])
#     first_bump_location = 0.5*(avg_bump_location_L[1] + avg_bump_location_R[1])
#     bump_locations[i] = last_bump_location - first_bump_location
# end 

# p = plot()
# Utils.plot_segments(J_head_range, bump_locations, tol = 3, xlabel = L"J", ylabel = L"\theta_{\text{bump mean}}")
# yticks!(
#     range(start = minimum(bump_locations), stop = maximum(bump_locations) , length = 5), 
#     [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]
# )
# vline!([-1, 1], legend = false, color = :black, linestyle = :dash)
# xticks!([-4, -2, -1, 0, 1, 2, 4], ["-4", "-2", "-1", "0", "1", "2", "4"])
# plot!(p, grid=true)
# sigma_str = replace(string(Q3.volatility), "." => "_")
# savefig("data/Q33_sigma_$sigma_str.pdf")

# # plot([-1.7, -1.7, 1.7, 1.7], [0, 0.5, 0.5, 0], color=:red)
# display(p)


# RandomTrajectory.plot_trajectory(pos[:,1], pos[:,2])