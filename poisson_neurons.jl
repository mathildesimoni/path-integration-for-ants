using LinearAlgebra

function sigmoid(z::Real)
    return 1.0 / (1.0 + exp(-z))
end

function g(x::Real, alpha::Real=1.0, beta::Real=0.0)
    return sigmoid(2.0 * alpha * (x - beta))
end

function h_for_loop(h_init::Real, delta_t::Real, n::Real, R::Real, I_t::Array, tau::Real)
    # I_t is the input already computed for every timestep (because otherwise we would have to pass the arguments the function I needs to this function too)
    h_t = zeros(Float64, n+1)
    h_t[1] = h_init
    for i in range(1, length=n)
        h_t[i+1] = (1 - delta_t/tau) * h_t[i] + (delta_t/tau) * R * I_t[i]
    end
    return h_t
end

function h(h_init::Real, delta_t::Real, n::Real, R::Real, I_t::Array, tau::Real)
    lower_diag = (delta_t/tau - 1) * ones(n)
    diag = ones(n+1)
    A = Bidiagonal(diag, lower_diag, :L)
    b = vcat([h_init], (delta_t * R / tau) * I_t[1:n])
    return A\b
end

function r(h::Real, alpha::Real=1.0, beta::Real=0.0, ro::Real=1.0)
    return ro * g(h, alpha, beta)
end

# slow oscillating input for the poisson neurons
function I(t::Real, Io::Real, omega::Real)
    return Io * sin(omega * t)
end

# Theoretical instantaneous firing rate
function r(t::Real, alpha::Real=1.0, beta::Real=0.0, ro::Real=1.0, Io::Real=1.0, omega::Real=1.0)
    return ro * g(I(t, Io, omega), alpha, beta)
end