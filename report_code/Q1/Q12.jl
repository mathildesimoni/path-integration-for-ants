using Revise
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons, SingleBumpAttractor, BumpAttractorUtils
using Random, Distributions
using Q1
using Utils



N = np.N
delta_t = sp.delta_t

n = sp.n
T = sp.T

x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)
h_init = rand(Uniform(0,1), N) # initial potential values sampled from the uniform distribution

# plot parameters
nb_ticks_x = 5

I_ext(x::Real, t::Real) = 0.0

spikes = SingleBumpAttractor.simulate_network(h_init, x_i, I_ext, 0.0, sp, np)
p = Utils.raster_plot(spikes, sp, np)
bin_size = 10 # ms, as advised in the instructions
Utils.plot_avg_bump_location(spikes, x_i, bin_size, sp, np, label="Average Bump Location", color=:orange)
savefig(p, "data/Q12_bump_location.pdf")

tau_values = [1.0, 5.0, 10.0]

for tau in tau_values
        np.tau = tau
        spikes = SingleBumpAttractor.simulate_network(h_init, x_i, I_ext, 0.0, sp, np)
        p = Utils.raster_plot(spikes, sp, np, title=LaTeXString("\$\\tau = $tau\$"))
        bin_size = 10 # ms, as advised in the instructions
        Utils.plot_avg_bump_location(spikes, x_i, bin_size, sp, np, label="Average Bump Location", color=:orange)
        tau_str = replace(string(tau), "." => "_")
        filename = "data/Q12_drift_tau_$(tau_str).pdf"
        savefig(p, filename)
end

