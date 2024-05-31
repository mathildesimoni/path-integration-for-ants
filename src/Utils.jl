module Utils
    using Statistics
    using BumpAttractorUtils
    using Plots, LaTeXStrings

    function split_into_segments(data::Array)
        diff_data = diff(data)
        diff_data = abs.(diff_data)
        tol = mean(diff_data) * 10
        idx = findall(x -> x > tol, diff_data)
        
        if length(idx) == 0
            return [data]
        end
        
        segments = []
        last_index = 1
        for i in eachindex(idx)
            push!(segments, data[last_index:idx[i]])
            last_index = idx[i] + 1
        end
        push!(segments, data[last_index:end])
        return segments
    end 

    function raster_plot(spikes::Array, sp::SimulationParameters, np::NetworkParameters)
        heatmap(transpose(spikes), 
                title="Network Activity", 
                xlabel=L"t"*" (ms)", 
                ylabel= "Neuron Location", 
                c = reverse(cgrad(:grayC)), 
                colorbar=false, 
                right_margin = 3Plots.mm, 
                left_margin = 2Plots.mm, 
                yticks = (range(start = 0, stop = np.N , length =5), 
                [L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2 \pi"]), 
                xticks = (Int.(0:sp.n/5:sp.n), 
                Int.(0:sp.T/5:sp.T)))
    end
end