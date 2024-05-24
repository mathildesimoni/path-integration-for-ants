module BumpAttractorUtils

    export locate_bump, locate_bump_avg, weighted_circular_mean, SimulationParameters, NetworkParameters

    using Parameters

    @with_kw mutable struct SimulationParameters
        T::Real=1000 
        n::Int=10_000
        delta_t::Real=0.1
        rate_neurons::Bool=false
        delay::Int=0
    end
    
    @with_kw mutable struct NetworkParameters
        N::Int=100
        R::Real=1.0
        tau::Real=10.0
        alpha::Real=2.0
        beta::Real=0.5
        ro::Real=1.0
        J::Real=1.0
    end
    
    function weighted_circular_mean(w, x)
        # implement the circular mean technique
        # https://en.wikipedia.org/wiki/Circular_mean#Using_complex_arithmetic
        cos_weighted_sum = sum(w .* cos.(x))
        sin_weighted_sum = sum(w .* sin.(x))
        return atan(sin_weighted_sum, cos_weighted_sum) 
    end

    function locate_bump(S_i::AbstractArray, x_i::Array)
        return mod(weighted_circular_mean(S_i, x_i), 2*pi) # within [0, 2pi] range
    end

    # the only difference with the previous function is the signature (need it to be different for broadcasting purposes)
    function locate_bump_avg(S_i::Array, x_i::AbstractArray)
        return mod(weighted_circular_mean(S_i, x_i), 2*pi) # within [0, 2pi] range
    end


end