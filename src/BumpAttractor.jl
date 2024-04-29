module BumpAttractor

    export m_cos, m_sin, I
    
    # first collective variable
    function m_cos(N::Int=100, x_i::Array{Float64}, S_i::Array{Int})
        return 1/N * sum(cos.(x_i) .* S_i)
    end

    # second collective variable
    function m_sin(N::Int=100, x_i::Array{Float64}, S_i::Array{Int})
        return 1/N * sum(sin.(x_i) .* S_i)
    end

    # Input function
    function I(J::Real, x::Real, m_cos::Real, m_sin::Real, I_ext::Real)
        return J * (cos(x) * m_cos + sin(x) * m_sin) + I_ext
    end

end