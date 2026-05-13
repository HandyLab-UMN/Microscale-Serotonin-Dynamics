# Consider a compartmental-reaction diffusion system with subsets of 
# compartments of radius εᵣ centred at {xj} in ℝ² with only one diffusing species u and intracellular species μ:
# extracellular diffusion & degradation
# ∂ₜ u = D Δu  - σ u,           t ∈ (0, ∞),    x ∈ ℝ²\⋃_{j=1}^n Bεᵣ(xj)
# reactive boundary conditions from intra- to/from extracellular molecules
# ε D ∂nj u = d₁ u - d₂ μj,    t ∈ (0, ∞),    x ∈ Bεᵣ(xj),    j ∈ {1, ..., n}
# intracellular reaction kinetics
# dₜ μj = f(μj) + ∑ₗ b*δ(t-Tlj) + 1/εᵣ ∫∂Bεᵣ(xj) (d₁ u - d₂ μj) dS
# Each compartment is located on an invisible 1-D neuronal fiber. If the corresponding fiber fires at Tlk, there is
# serotonin release into all compartments along the fiber with the kth compartment. The index l of Tlk ranges over
# all firing instances of the fiber with the kth compartment. The Tlj are random times at which the concentration 
# in the jth compartment increases by b.
#
# This code numerically solves the corresponding asymptotically reduced integro-ODE system
# dₜ μj = f(μj) + ∑ₗ c*δ(t-Tlj) + Bj(t)
# ∫₀ᵗ Bj'(τ) E₁(σ(t-τ)) dτ = αj Bj + γj μj + ∑{k≠j} ∫₀ᵗ Bk(τ)/(t-τ) exp(-|xj-xk|²/(4D(t-τ)) - σ(t-τ)) dτ
# for j ∈ {1, ..., n}
#
# using an IMEX method for implicitly (with Backward-Euler method) integrating
# dₜ μjB = Bj(t) and explicitly (with RK4) integrating the rest (except for Bj),
#  i.e., dₜ μj = f(μj, ηj) ∧ dₜ ηj = g(μj, ηj), using the derived scheme for Bj.

# Units:
# 1 unit in x-direction means 5 μm 
# 1 unit in t-direction means 1/62.5 s
# D = 1 means 1000 μm²/s
# uj(t) = 1 means 10^{-22} mol 
# U(t, x) = 1 means 1/25*10^{-22} mol/μm²


include("PackagesAndParameters.jl")

tarray = readdlm("./datasets/tarray60s.txt")
μarray = readdlm("./datasets/mu_array_freq$(freq)_$(n_on_fiber)line_$(modelType).txt")
μarray = reshape(μarray, (length(tarray), ncells))


indcs = 5000:length(tarray)  # 179660 # 6271 # 179660:179660+6271
colors = [:magenta, :dodgerblue, :indianred, :orange, :purple, :orchid, :blue, :red, :darkorange1, :rebeccapurple]

# spike maximum plot
maxμarray = zeros(ntpoints, ncells)
for tidx in 1:ntpoints
	if mod(tidx,1000)==0
    	println("max t-loop: \t $(round(Int, tidx/ntpoints * 100))%")
	end
   maxμarray[tidx, :] .= maximum(μarray[1:tidx, :], dims=1)'
end

μeq = steadystate(freq_nondim, 4. .* sum(firingarray[fireidx, :] for fireidx in 1:nfiringfibers) .+ 1,modelType)  # [2., 1.5, 1., 1.5, 2.]

maxapprox = μeq .+ get_spikeheight(μeq, std_kick)

writedlm("./datasets/max_mu$(freq)_$(n_on_fiber)line.txt", maxμarray)
writedlm("./datasets/max_approx$(freq)_$(n_on_fiber)line.txt", maxapprox)
writedlm("./datasets/tarray$(freq)_maxCompare.txt", tarray./kR)
writedlm("./datasets/indcs$(freq)_maxCompare.txt", indcs)
