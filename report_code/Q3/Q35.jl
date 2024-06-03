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

Q3.init_params()

np.tau = 10.0
sp.delta_t = 0.1
sp.T = 3000
sp.n = Int64(sp.T/sp.delta_t)
np.N = 1000

angles, pos = RandomTrajectory.random_trajectory(Q3.speed, sp.T, sp.delta_t, Q3.volatility)

x_i = collect(range(start = 0, stop = 2*pi, length = np.N + 1)[1:np.N]) # equally spaced neurons over the range [0, 2pi)
h_init = rand(Uniform(0,1), np.N) # initial potential values sampled from the uniform distribution

sp.delay = 500

t_to_idx(t) = Int(floor(t/sp.delta_t)) + 1
theta(t) = angles[t_to_idx(t)]
I_ext(x, t) = Q3.I_ext_head(x, t, Q3.Io, theta)

np.J = 5
# sp.rate_neurons = true
spikes = SingleBumpAttractor.simulate_network(h_init, x_i, I_ext, 0.0, sp, np)

# J_head = 0.85
# J_head = 0.35
# J_head = 1.5
J_head = 0.5

I_head_val_cos = (spikes/sp.delta_t) * cos.(x_i)
I_head_val_sin = (spikes/sp.delta_t) * sin.(x_i)
I_head_x(t) = I_head_val_cos[Q3.t_to_idx(t)]
I_head_y(t) = I_head_val_sin[Q3.t_to_idx(t)]

np.J = 3
spikes_L_x, spikes_R_x = Q3.integrate_position(J_head, I_head_x)
spikes_L_y, spikes_R_y = Q3.integrate_position(J_head, I_head_y)

spikes_x = (spikes_L_x + spikes_R_x)./2
spikes_y = (spikes_L_y + spikes_R_y)./2

bin_size = 10
bump_location_x = Utils.spikes_to_average_bump_location(spikes_x, x_i, bin_size, sp)
bump_location_y = Utils.spikes_to_average_bump_location(spikes_y, x_i, bin_size, sp)

# Plot trajectory
plot()
RandomTrajectory.plot_trajectory(pos[:,1], pos[:,2], label = "Trajectory", start_end_labels = true, color=:blue)
savefig("data/Q35_trajectory.svg")

# Raw results
function plot_raw_results()
    p = plot()
    RandomTrajectory.plot_trajectory(pos[:,1], pos[:,2], label = "Trajectory", start_end_labels = false, color=:blue)
    RandomTrajectory.plot_trajectory(bump_location_x, bump_location_y, label = "Integrated", start_end_labels = false, color=:black, alpha=0.15)
    smoothed_bump_location_x = Utils.smooth(bump_location_x, 10)
    smoothed_bump_location_y = Utils.smooth(bump_location_y, 10)
    RandomTrajectory.plot_trajectory(smoothed_bump_location_x, smoothed_bump_location_y, label = "Smoothed", start_end_labels = true, color=:orange)

    
    # vline and hline and x,y = pi
    vline!([pi], color = :black, lw = 1, label = false, linestyle = :dash)
    hline!([pi], color = :black, lw = 1, label = false, linestyle = :dash)
    xticks!([0, 1, 2, pi], ["0", "1", "2", L"\pi"])
    yticks!([0, 1, 2, pi], ["0", "1", "2", L"\pi"])
    savefig("data/Q35_raw_results.svg")
    return p
end
# display(plot_raw_results())

# Manipulate results
bump_location_x = bump_location_x .- bump_location_x[1]
bump_location_y = bump_location_y .- bump_location_y[1]

# Plot x 
function plot_x_results()
    # raster plot
    p = Utils.raster_plot(spikes_x, sp, np)
    bin_size = 10
    Utils.plot_avg_bump_location(spikes_x, x_i, bin_size, sp, np, color=:orange, label="Integrated")

    target = pos[:, 1]
    target .+= pi
    target = target .% (2*pi)
    # plot!(Utils.map_angle_to_idx.(target, np.N), label="Target", color=:blue)
    ylabel!("x")
    xlabel!("Time step")
    savefig("data/Q35_x_results.svg")
    return p
end
# display(plot_x_results())

function plot_x_results_no_raster()
    plot(pos[:,1], label="Target", color=:blue)
    bin_length = Int64(bin_size/sp.delta_t)
    ylabel!("x")
    xlabel!("Time step")
    return plot!(0:bin_length:sp.n-1, bump_location_x, label="Integrated", color=:orange)
end
# display(plot_x_results_no_raster())

function plot_y_results_no_raster()
    plot(pos[:,2], label="Target", color=:blue)
    bin_length = Int64(bin_size/sp.delta_t)
    ylabel!("y")
    xlabel!("Time step")
    return plot!(0:bin_length:sp.n-1, bump_location_y, label="Integrated", color=:orange)
end
# display(plot_y_results_no_raster())

function plot_y_results()
    # raster plot
    p = Utils.raster_plot(spikes_y, sp, np)
    bin_size = 10
    Utils.plot_avg_bump_location(spikes_y, x_i, bin_size, sp, np, color=:orange, label="Integrated")

    target = pos[:, 2]
    target .+= pi
    target = target .% (2*pi)
    # plot!(Utils.map_angle_to_idx.(target, np.N), label="Target", color=:blue)
    ylabel!("y")
    xlabel!("Time step")
    savefig("data/Q35_y_results.svg")
    return p
end

display(plot_y_results())

function plot_centered_results()
    p = plot()
    RandomTrajectory.plot_trajectory(pos[:,1], pos[:,2], label = "Trajectory", start_end_labels = false, color=:blue)
    RandomTrajectory.plot_trajectory(bump_location_x, bump_location_y, label = "Integrated", start_end_labels = false, color=:black, alpha=0.15)

    # smoothed version
    window_size = 10
    smoothed_bump_location_x = Utils.smooth(bump_location_x, window_size)
    smoothed_bump_location_y = Utils.smooth(bump_location_y, window_size)
    RandomTrajectory.plot_trajectory(smoothed_bump_location_x, smoothed_bump_location_y, label = "Smoothed", start_end_labels = true, color=:orange)
    savefig("data/Q35_results_scaled.svg")
    return p
end

function plot_scaled_results()
    p = plot()
    RandomTrajectory.plot_trajectory(pos[:,1], pos[:,2], label = "Trajectory", start_end_labels = false, color=:blue)
    
    scaled_bump_location_x = bump_location_x .* 0.55
    # scaled_bump_location_x = bump_location_x
    scaled_bump_location_y = bump_location_y .* 0.55
    # scaled_bump_location_y = bump_location_y 
    

    RandomTrajectory.plot_trajectory(scaled_bump_location_x, scaled_bump_location_y, label = "Integrated", start_end_labels = false, color=:black, alpha=0.15)
    smoothed_bump_location_x = Utils.smooth(scaled_bump_location_x, 5)
    smoothed_bump_location_y = Utils.smooth(scaled_bump_location_y, 5)
    RandomTrajectory.plot_trajectory(smoothed_bump_location_x, smoothed_bump_location_y, label = "Smoothed", start_end_labels = true, color=:orange)
    savefig("data/Q35_volatile_results_scaled.svg")
    return p
end
display(plot_scaled_results())


