using Revise
a = push!(LOAD_PATH, pwd()*"/src")
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons, SingleBumpAttractor
using Random, Distributions
using Q1

# default parameters
N = np.N
n = sp.n
T = sp.T
delta_t = sp.delta_t

x_i = collect(range(start = 0, stop = 2*pi, length = N + 1)[1:N]) # equally spaced neurons over the range [0, 2pi)
h_init = rand(Uniform(0,1), N) # initial potential values sampled from the uniform distribution
I_ext_bool = false # no external input
I_ext(x::Real, t::Real) = 0.0

# plot parameters
nb_ticks_x = 5
bin_size = 10 # ms, as advised in the instructions

np_test = deepcopy(np)
sp_test = deepcopy(sp)

# Ranging number of poissons neurons
N_values = [50, 100, 300, 500, 1000]
for N_test in N_values

    np_test.N = N_test
    x_i_test = collect(range(start = 0, stop = 2*pi, length = N_test + 1)[1:N_test])
    h_init_test = rand(Uniform(0,1), N_test)

    # simulate the network activity
    spikes = SingleBumpAttractor.simulate_network(h_init_test, x_i_test, I_ext, 0.0, sp, np_test)
    
    # find the location of the bump at every timestep and average over small bins of time
    bump_location = locate_bump.(eachrow(spikes), Ref(x_i_test))
    bin_length = Int64(bin_size/delta_t)
    bump_location_bins = transpose(reshape(bump_location[1:n], bin_length, Int((n)/bin_length)))
    avg_bump_location = locate_bump_avg.(Ref(ones(bin_length)), eachrow(bump_location_bins)) # need to use a circular mean method again
    heatmap(transpose(spikes), 
            title="Network Activity (" * string(N_test) * " Neurons)", 
            xlabel=L"t"*" (ms)", 
            ylabel= "Neuron Location", 
            c = reverse(cgrad(:grayC)), 
            colorbar=false, 
            right_margin = 3Plots.mm, 
            left_margin = 2Plots.mm, 
            yticks = (range(start = 0, stop = N_test , length =5), 
            [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), 
            xticks = (Int.(0:n/nb_ticks_x:n), 
            Int.(0:T/nb_ticks_x:T))
            )
    plot!(0:bin_length:n-1, avg_bump_location * (N_test/(2*pi)), label = "center of the bump")
    filename = "data/Q13_N_" * string(N_test) * ".pdf"
    savefig(filename)
end

np_test.N = np.N # reset to default value

# Ranging timescale
tau_values = [1, 5, 10, 20, 30, 40]
for tau_test in tau_values
    
    np_test.tau = tau_test

    # simulate the network activity
    spikes = SingleBumpAttractor.simulate_network(h_init, x_i, I_ext, 0.0, sp, np_test)

    # find the location of the bump at every timestep and average over small bins of time
    bump_location = locate_bump.(eachrow(spikes), Ref(x_i))
    bin_length = Int64(bin_size/delta_t)
    bump_location_bins = transpose(reshape(bump_location[1:n], bin_length, Int((n)/bin_length)))
    avg_bump_location = locate_bump_avg.(Ref(ones(bin_length)), eachrow(bump_location_bins)) # need to use a circular mean method again

    heatmap(transpose(spikes), 
            title="Network Activity (τ = " * string(tau_test) * ")", 
            ylabel= "Neuron Location", 
            c = reverse(cgrad(:grayC)), 
            colorbar=false, 
            right_margin = 3Plots.mm, 
            left_margin = 2Plots.mm, 
            yticks = (range(start = 0, stop = N , length =5), 
            [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), 
            xticks = (Int.(0:n/nb_ticks_x:n), 
            Int.(0:T/nb_ticks_x:T))
            )
    plot!(0:bin_length:n-1, avg_bump_location * (N/(2*pi)), label = "center of the bump")
    filename = "data/Q13_tau_" * string(tau_test) * ".pdf"
    savefig(filename)
end

np_test.tau = np.tau # reset to default value


# Ranging timestep
delta_t_values = [0.05, 0.1, 0.2, 0.25, 0.5]
for delta_t_test in delta_t_values

    sp_test.delta_t = delta_t_test
    n_test = Int64(T/delta_t_test)
    sp_test.n = n_test

    # simulate the network activity
    spikes = SingleBumpAttractor.simulate_network(h_init, x_i, I_ext, 0.0, sp_test, np)

    # find the location of the bump at every timestep and average over small bins of time
    bump_location = locate_bump.(eachrow(spikes), Ref(x_i))
    bin_length = Int64(bin_size/delta_t_test)
    bump_location_bins = transpose(reshape(bump_location[1:n_test], bin_length, Int((n_test)/bin_length)))
    avg_bump_location = locate_bump_avg.(Ref(ones(bin_length)), eachrow(bump_location_bins)) # need to use a circular mean method again

    heatmap(transpose(spikes), 
            title="Network Activity (Δ = " * string(delta_t_test) * ")", 
            xlabel=L"t"*" (ms)", 
            ylabel= "Neuron Location", 
            c = reverse(cgrad(:grayC)), 
            colorbar=false, 
            right_margin = 3Plots.mm, 
            left_margin = 2Plots.mm, yticks = (range(start = 0, stop = N , length =5), 
            [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), 
            xticks = (Int.(0:n_test/nb_ticks_x:n_test), 
            Int.(0:T/nb_ticks_x:T))
            )
    plot!(0:bin_length:n_test-1, avg_bump_location * (N/(2*pi)), label = "center of the bump")
    filename = "data/Q13_delta_t_" * string(delta_t_test) * ".pdf"
    savefig(filename)
end

sp_test.delta_t = sp.delta_t
sp_test.n = sp.n





