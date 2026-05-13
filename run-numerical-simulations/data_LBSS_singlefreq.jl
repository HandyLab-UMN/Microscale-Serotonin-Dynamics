# Finds the periodic steady state and the rolling average 
# Must run run_simulation_2D first!

include("PackagesAndParameters.jl")

tarray = readdlm("./datasets/tarray60s.txt")
μarray = readdlm("./datasets/mu_array_freq$(freq)_$(n_on_fiber)line_$(modelType).txt")
μarray = reshape(μarray, (length(tarray), ncells))

# adding steady-state to the plot # 614500
μeq = steadystate(freq_nondim, 4. .* sum(firingarray[fireidx, :] for fireidx in 1:nfiringfibers) .+ 1,modelType)  # [2., 1.5, 1., 1.5, 2.]
writedlm("./datasets/mu_eq$(freq)_$(n_on_fiber)line_$(modelType).txt", μeq)

# adding the true averaged dynamics
ntsteps_period = Int(floor(1/freq_nondim / Δt_finer))
tidx_end_ave = ntpoints-ntsteps_period
μave_array = zeros(tidx_end_ave, 5)
for j in 1:5
    for tidx in 1:tidx_end_ave
        μave_array[tidx, j] = mean(μarray[tidx:tidx+ntsteps_period, j])
    end
end
writedlm("./datasets/mu_ave$(freq)_$(n_on_fiber)line_$(modelType).txt", μave_array)





