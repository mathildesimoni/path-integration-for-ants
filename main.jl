# will take a little moment to install the packages on the first run
using Pkg
Pkg.add("Plots")
Pkg.add("LaTeXStrings") 

using Plots, LaTeXStrings

# include the functions from the other files
include("poisson_neurons.jl")

# define the parameters
h = range(-5, 5, length=1000)
alpha_values = [0.5 1.0 2.0 100.0]
beta_values = [-4.0 0.0 4.0]


# ranging alpha
beta = 0.0
p1 =plot()
# plot!(title = "Activation \$g\$ vs. Potential \$h\$", xlabel=L"x", ylabel=L"g(x)")
plot!(xlabel=L"h", ylabel=L"g")
for alpha in alpha_values
    # the dot . here is used to broadcast the function g over the array h (g is defined for scalars)
    result = g.(h, alpha, beta)
    plot!(h, result, label=L"\alpha=%$(alpha)", lw=2)
end

# ranging beta
p2 =plot()
alph = 1.0
plot!(xlabel=L"h", ylabel=L"g")
for beta in beta_values
    result = g.(h, alph, beta)
    plot!(h, result, label=L"\beta=%$(beta)", lw=2)
end

# plot(p1, p2, layout=(2,1), size=(800, 800))
plot(p1)
# plot(p2)


# save
savefig("data/sigmoid_alpha.png")
# savefig("data/sigmoid_beta.png")

# plot(h, g_values, label=labels, lw=2)

