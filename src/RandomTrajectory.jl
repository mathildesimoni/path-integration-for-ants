module RandomTrajectory

    using Distributions

    function random_trajectory(
                speed::Real=1.0, # speed of the ant (m/s) (if 1, then the ant moves 1m in 1000ms)
                T::Real=1000,  # total time of trajectory (ms)
                delta_t::Real=0.1, # time step (ms)
                volatility::Real=1.0,
    )
        n = Int64(round(T/delta_t))
        pos = zeros(Float64, (n+1, 2))
        angles = zeros(Float64, n+1)
        angles[1] = pi
        for i in range(1, length=n)
            angles[i+1] = rand(Normal(angles[i], volatility))
            pos[i+1, :] = pos[i, :] + speed *  delta_t * [cos(angles[i+1]), sin(angles[i+1])]
        end
        return angles, pos
    end

    export random_trajectory

end 