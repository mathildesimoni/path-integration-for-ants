using Revise
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Distributions
println(Revise.errors())

using BumpAttractorUtils, SingleBumpAttractor
using Q3
using RandomTrajectory
using Utils

N = np.N
n = sp.n
T = sp.T
delta_t = sp.delta_t
bin_size = 5 # ms, as advised in the instructions
bin_length = Int64(bin_size/sp.delta_t)
x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)
h_init = rand(Uniform(0,1), N) # initial potential values sampled from the uniform distribution
sp.delay = 500

Io = 1
sigmas = [0.01, 0.08]

for sigma in sigmas
    
    # external input
    angles, pos = RandomTrajectory.random_trajectory(Q3.speed, T, delta_t, sigma)
    plot(title = "Ant Trajectory")
    RandomTrajectory.plot_trajectory(pos[:,1], pos[:,2])
    savefig("data/Q32_trajectory_sigma_" * string(sigma) * ".pdf")

    theta(t) = angles[Q3.t_to_idx(t)]
    I_ext(x, t) = Q3.I_ext_head(x, t, Io, theta)

    # simulate network
    S_i = SingleBumpAttractor.simulate_network(h_init, x_i, I_ext, 0.0, sp, np)

    # plot
    p = Utils.raster_plot(S_i, sp, np, title = "Bump Location and Head Direction")
    Utils.plot_angle_location(angles, N, color = :lightblue, label=LaTeXString(L"\theta_{\mathrm{input}}^{\mathrm{H}}"), tol = 30)
    Utils.plot_avg_bump_location(S_i, x_i, bin_size, sp, np, color = :red, label=LaTeXString(L"\theta_{\mathrm{bump}}^{\mathrm{H}}"), tol = 10)
    display(p)
    savefig("data/Q32_bump_sigma_" * string(sigma) * ".pdf")
end