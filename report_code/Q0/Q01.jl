a = push!(LOAD_PATH, pwd()*"/src")
using Plots, LaTeXStrings
default(fontfamily="Computer Modern")
using Neurons

# define the parameters
h_array = range(-5, 5, length=1000)
alpha_values = [0.5 1.0 2.0 100.0]
beta_values = [-2.0 0.0 2.0]

# ranging alpha
beta = 0.0
p1 = plot()
# plot!(title = "Activation \$g\$ vs. Potential \$h\$", xlabel=L"x", ylabel=L"g(x)")
plot!(xlabel=L"h", ylabel=L"g")
for alpha in alpha_values
    # the dot . here is used to broadcast the function g over the array h_arr (g is defined for scalars)
    result = g.(h_array, alpha, beta)
    plot!(h_array, result, label=L"\alpha=%$(alpha)", lw=2)
end

# ranging beta
p2 =plot()
alph = 1.0
plot!(xlabel=L"h", ylabel=L"g")
for beta in beta_values
    result = g.(h_array, alph, beta)
    plot!(h_array, result, label=L"\beta=%$(beta)", lw=2)
end

plot(p1)
savefig("data/Q01_sigmoid_alpha.pdf")
plot(p2)
savefig("data/Q01_sigmoid_beta.pdf")
