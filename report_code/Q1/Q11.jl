using Revise
a = push!(LOAD_PATH, pwd()*"/src", @__DIR__)
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons, SingleBumpAttractor, BumpAttractorUtils
using Random, Distributions
using Q1

N = np.N
J = np.J
n = sp.n
T = sp.T

x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)
h_init = rand(Uniform(0,1), N) # initial potential values sampled from the uniform distribution

# plot parameters
nb_ticks_x = 5

I_ext(x::Real, t::Real) = 0.0

filename = "data/Q11_J_" * string(J) * ".pdf"
spikes = SingleBumpAttractor.simulate_network(h_init, x_i, I_ext, 0.0, sp, np)
heatmap(transpose(spikes), title="Network Activity", xlabel=L"t"*" (ms)", ylabel= "neuron location", c = :grayC, colorbar=false, right_margin = 5Plots.mm, left_margin = 2Plots.mm, yticks = (range(start = 0, stop = N , length =5), [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), xticks = (Int.(0:n/nb_ticks_x:n), Int.(0:T/nb_ticks_x:T)))

# ranging J
J_values = [0.1, 0.5, 1, 2, 3, 4, 5, 6, 10, 20, 50, 100]
for J in J_values
    np.J = J
    filename = "data/Q11_J_" * string(J) * ".pdf"
    spikes = SingleBumpAttractor.simulate_network(h_init, x_i, I_ext, 0.0, sp, np)
    heatmap(transpose(spikes), 
            title = "Network Activity", 
            xlabel = L"t"*" (ms)", 
            ylabel = "neuron location", 
            c = reverse(cgrad(:grayC)), 
            colorbar = false, 
            right_margin = 5Plots.mm, 
            left_margin = 2Plots.mm, 
            yticks = (range(start = 0, stop = np.N , length =5), 
            [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), 
            xticks = (Int.(0:n/nb_ticks_x:n), 
            Int.(0:T/nb_ticks_x:T)))
    savefig(filename) 
end