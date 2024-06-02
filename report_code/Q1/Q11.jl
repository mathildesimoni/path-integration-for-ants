using Revise
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons, SingleBumpAttractor, BumpAttractorUtils
using Random, Distributions
using Q1
using Utils

N = np.N
J = np.J
n = sp.n
T = sp.T

x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)
h_init = rand(Uniform(0,1), N) # initial potential values sampled from the uniform distribution

# plot parameters
nb_ticks_x = 5

I_ext(x::Real, t::Real) = 0.0

# ranging J
J_values = [0, 4.6, 4.7, 4.8, 5]
for J in J_values
    np.J = J
    str_J = replace(string(J), "." => "_")
    filename = "data/Q11_J_" * str_J * ".pdf"
    spikes = SingleBumpAttractor.simulate_network(h_init, x_i, I_ext, 0.0, sp, np)
    Utils.raster_plot(spikes, sp, np, title="")
end

np.J = 5
spikes = SingleBumpAttractor.simulate_network(h_init, x_i, I_ext, 0.0, sp, np)
str_J = replace(string(np.J), "." => "_")
Utils.raster_plot(spikes, sp, np, title=LaTeXString("\$J=$(np.J)\$"))
filename = "data/Q11_J_" * str_J * ".pdf"
savefig(filename) 