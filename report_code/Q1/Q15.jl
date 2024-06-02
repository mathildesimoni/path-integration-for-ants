using Revise
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons, SingleBumpAttractor
using Random, Distributions
using Q1
using Utils

N = np.N
n = sp.n
T = sp.T
delta_t = sp.delta_t

x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)
h_init = rand(Uniform(0,1), N) # initial potential values sampled from the uniform distribution

# plot parameters
nb_ticks_x = 5

# no external input
I_ext(x::Real, t::Real) = 0.0

# phi_values = [0.01, 0.05, -0.01]
phi_values = [-0.01]
for phi in phi_values
        # simulate the network activity
        spikes = SingleBumpAttractor.simulate_network(h_init, x_i, I_ext, phi, sp, np)

        # average location of the bump over small bins of time (to have a clearer plot)
        p = Utils.raster_plot(spikes, sp, np, title=LaTeXString("\$\\varphi=$phi\$"))
        bin_size = 10 # ms, as advised in the instructions
        Utils.plot_avg_bump_location(spikes, x_i, bin_size, sp, np, tol=5)
        phi_str = replace(string(phi), "." => "_")
        filename = "data/Q15_phi_" * phi_str * ".pdf"
        savefig(filename)
end
