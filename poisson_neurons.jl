function sigmoid(z::Real)
    return 1.0 / (1.0 + exp(-z))
end

function g(x::Real, alpha::Real=1.0, beta::Real=0.0)
    return sigmoid(2.0 * alpha * x - beta)
end

