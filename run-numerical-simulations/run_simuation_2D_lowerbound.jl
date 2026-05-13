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

# time domain
tfinal_extended = tfinal+800.  
tarray = vcat([0, Δt], collect(2*Δt:Δt_finer:tfinal_extended))
ntpoints = length(tarray)  # == tfinal/Δt + 1

Tfiber = []
nfirings = round(Int, (tfinal_extended- 0.6*kR) * freq_nondim)  # no switching between frequencies
for fiberidx in firingfibers
	append!(Tfiber, [[0.6*kR + (k-1) * 1/freq_nondim for k in 1:nfirings]])  # no switching between frequencies
end


fiberindex = firingfibers[1]
periodarray = Tfiber[fiberindex][2:end] .- Tfiber[fiberindex][1:end-1]
kick_averaged_array = [1/P for P in periodarray]
kick_averaged_array = vcat([0, 0, 0], kick_averaged_array, [0, 0, 0, 0])
tkicks = vcat([j/5 * Tfiber[fiberindex][1] for j in 1:3], Tfiber[fiberindex], [Tfiber[fiberindex][end] + j/5 * (tfinal_extended - Tfiber[fiberindex][end]) for j in 3:5])
kick = Spline1D(tkicks, kick_averaged_array; k=1, s=0)  
function fvecLB(t::Float64, μ::Array{Float64, 1}) :: Array{Float64, 1}
    # on the fly time-averaged f for switching of firing frequencies
    return -V .* μ ./ (K .+ μ) .+ b * kick(t) .* sum(firingarray[fireidx, :] for fireidx in 1:nfiringfibers)
end


#setting up a matrix with all weight arrays for all cell posiitions
wmatrix = zeros(Complex, ncells, ncells, 2*m+1)
for i in 1:ncells
    println("Computing exponential sum weights for the heat kernel: \t $(round(Int, i/ncells * 100))%")
    for j in i+1:ncells
        wmatrix[i, j, :] .= warray(norm(cpos[i] .- cpos[j], 2))
        wmatrix[j, i, :] .= wmatrix[i, j, :]
    end
end

M = m
@assert M == m "A new sarray of length 2*M+1 has to be created."
earray = [λ*h/(2*π) * cos(α + complex(0,1)*k*h) * log(1 + sarray[k+M+1]/σ)/sarray[k+M+1] for k in -M:M]

# creating the matrices for the time-stepping
b1 = Δt/2 * expint(1, σ*Δt) + (1-exp(-σ*Δt))/σ + 1/(2*σ)*exp(-σ*Δt) - 1/(2*σ^2*Δt) * (1-exp(-σ*Δt))
b2 = Δt/2 * expint(1, σ*Δt) - 1/(2*σ)*exp(-σ*Δt) + 1/(2*σ^2*Δt) * (1-exp(-σ*Δt))
b3 = zeros(Complex, length(sarray))
b4 = zeros(Complex, length(sarray))
b40 = zeros(Complex, length(sarray))  # for second time step
for l in eachindex(sarray)
    sl = sarray[l]
    b3[l] = exp(sl*Δt) * (-1/sl + 1/sl^2 * (exp(sl*Δt) - 1)/Δt)
    b4[l] = exp(sl*Δt) * ((exp(sl*Δt) - 1)/sl + 1/sl - 1/sl^2 * (exp(sl*Δt) - 1) / Δt)
    b40[l] = exp(sl*Δt) * ((exp(sl*Δt) - 1)/sl)
end

# with smaller time steps after the first two
b1_finer = Δt_finer/2 * expint(1, σ*Δt_finer) + (1-exp(-σ*Δt_finer))/σ + 1/(2*σ)*exp(-σ*Δt_finer) - 1/(2*σ^2*Δt_finer) * (1-exp(-σ*Δt_finer))
b2_finer = Δt_finer/2 * expint(1, σ*Δt_finer) - 1/(2*σ)*exp(-σ*Δt_finer) + 1/(2*σ^2*Δt_finer) * (1-exp(-σ*Δt_finer))
b3_finer = zeros(Complex, length(sarray))
b4_finer = zeros(Complex, length(sarray))
for l in eachindex(sarray)
    sl = sarray[l]
    b3_finer[l] = exp(sl*Δt_finer) * (-1/sl + 1/sl^2 * (exp(sl*Δt_finer) - 1)/Δt_finer)
    b4_finer[l] = exp(sl*Δt_finer) * ((exp(sl*Δt_finer) - 1)/sl + 1/sl - 1/sl^2 * (exp(sl*Δt_finer) - 1) / Δt_finer)
end


# solution arrays
μarray = zeros(ntpoints, ncells)
Barray = zeros(ntpoints, ncells)

# updating matrices for B
A = zeros(ncells, ncells)
A .= Diagonal(b1 .- αvec .* Δt)
M1 = zeros(ncells, ncells)
m1diag = (b1 - Δt * sum(earray .* exp.(2*Δt .* sarray)))
@assert imag(m1diag) < 1e-14
M1 .= Diagonal(real(m1diag) .* ones(ncells))

for j in 1:ncells
    for k in 1:ncells
        if j != k
            absxsq = ((cpos[j][1] - cpos[k][1])^2 + (cpos[j][2] - cpos[k][2])^2) / D
            A[j, k] = -Δt * expint(1, absxsq/(4*Δt))  # correction terms as well?
            m1 =  Δt * sum(wmatrix[j, k, :] .* b40)
            @assert abs(imag(m1)) < 1e-8 "abs(imag(m1)) = $(abs(imag(m1)))"
            M1[j, k] = real(m1)
        end
    end
end

# updating matrices for B after the first two time steps
A_finer = zeros(ncells, ncells)
A_finer .= Diagonal(b1_finer .- αvec .* Δt_finer)
ℳ_finer = zeros(ncells, ncells)
ℳdiag_finer = (b1_finer - b2_finer - sum(earray .* b3_finer)).*ones(ncells)
@assert maximum(abs.(imag.(ℳdiag_finer))) < 1e-8 "maximum(abs.(imag.(ℳdiag_finer))) = $(maximum(abs.(imag.(ℳdiag_finer))))"
ℳ_finer .= Diagonal(real.(ℳdiag_finer))
𝒩_finer = zeros(ncells, ncells)
𝒩diag_finer = (b2_finer + sum(earray .* (b3_finer .- b4_finer))).*ones(ncells)
@assert maximum(abs.(imag.(𝒩diag_finer))) < 1e-8 "maximum(abs.(imag.(𝒩diag_finer))) = $(maximum(abs.(imag.(𝒩diag_finer))))"
𝒩_finer .= Diagonal(real.(𝒩diag_finer))
𝒩_finer3 = zeros(ncells, ncells)
𝒩diag_finer3 = (b2_finer + sum(earray .* b3_finer)).*ones(ncells)
@assert maximum(abs.(imag.(𝒩diag_finer3))) < 1e-8 "maximum(abs.(imag.(𝒩diag_finer3))) = $(maximum(abs.(imag.(𝒩diag_finer3))))"
𝒩_finer3 .= Diagonal(real.(𝒩diag_finer3))
for j in 1:ncells
    for k in 1:ncells
        if j != k
            absxsq = ((cpos[j][1] - cpos[k][1])^2 + (cpos[j][2] - cpos[k][2])^2) / D
            A_finer[j, k] = -Δt_finer * expint(1, absxsq/(4*Δt_finer))  # correction terms as well?
            m1 =  Δt_finer * sum(wmatrix[j, k, :] .* b3_finer)
            @assert abs(imag(m1)) < 1e-8 "abs(imag(m1)) = $(abs(imag(m1)))"
            ℳ_finer[j, k] = real(m1)
            nval =  Δt_finer * sum(wmatrix[j, k, :] .* b4_finer)
            @assert abs(imag(nval)) < 1e-8 "abs(imag(nval)) = $(abs(imag(nval)))"
            𝒩_finer[j, k] = real(nval)
            𝒩_finer3[j, k] = real(nval)
        end
    end
end

# initial concentration in all cells, no firing yet
μarray[1, :] .= 0. * ones(ncells) # .+ pert
Barray[1, :] .= zeros(ncells)
Barray[2, :] .= -4*π*D ./ log.(Δt./(εᵣ^2 * exp(-γe) .* κ0vec)) .* βvec .* μarray[1, :]
dB2 = 4*π*D ./ (Δt .* log.(Δt./(εᵣ^2 * exp(-γe) .* κ0vec)).^2) .* βvec .* μarray[1, :]

Hw = zeros(Complex, 2*m+1, ncells)
He = zeros(Complex, 2*M+1, ncells)
F1 = zeros(Complex, ncells)
F2 = zeros(Complex, ncells)

### Instead with increased fineness
μarray2 = zeros(ncells)
μarray2 .= μarray[1, :]
for tidx1 in 1:tfineness_factor
    global μarray2
    # Rk4 weights
    t = (tidx1-1) * Δt_finer

    k1 = fvecLB(t, μarray2)
    k2 = fvecLB(t + Δt/2, μarray2 .+ Δt_finer/2 .* k1)
    k3 = fvecLB(t + Δt/2, μarray2 .+ Δt_finer/2 .* k2)
    k4 = fvecLB(t + Δt, μarray2 .+ Δt_finer .* k3)

    μarray2 .= μarray2 .+ Δt_finer/6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4) + Δt_finer .* Barray[2, :]  # SEMIEXPLICIT
end
μarray[2, :] += μarray2  # WARNING: "+=" is important due to the incorporation of fiber firing serotonin-release

### Instead with increased fineness
Barray[2, :] .= -4*π*D ./ log.(Δt./(εᵣ^2 * exp(-γe) .* κ0vec)) .* βvec .* μarray[2, :] .* (1 .- π^2 ./ (6 .* log.(Δt./(εᵣ^2 * exp(-γe) .* κ0vec)).^2))

μarray3 = zeros(ncells)
μarray3 .= μarray[2, :]
for tidx2 in 1:tfineness_factor
    global μarray3
    # Rk4 weights
    t = Δt + (tidx2-1) * Δt_finer
    
    k1 = fvecLB(t, μarray3)
    k2 = fvecLB(t + Δt/2, μarray3 .+ Δt_finer/2 .* k1)
    k3 = fvecLB(t + Δt/2, μarray3 .+ Δt_finer/2 .* k2)
    k4 = fvecLB(t + Δt, μarray3 .+ Δt_finer .* k3)

    μarray3 .= μarray3 .+ Δt_finer/6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4) + Δt_finer .* Barray[2, :]  # SEMIEXPLICIT
end
μarray[3, :] += μarray3
# SEMIEXPLICIT >>>
Bnext = A \ (M1*Barray[2, :] .+ Δt .* Γ*μarray[3, :] .- Δt * b2 .* dB2)
@assert maximum(abs.(imag.(Bnext))) < 1e-8 "maximum(abs.(imag.(Bnext))) = $(maximum(abs.(imag.(Bnext))))"
Barray[3, :] .= real.(Bnext)
# <<<

# updating the history terms for the second time step
if ncells > 1
    Hw .= b40 * transpose(Barray[2, :])
end
He .= exp.(2*Δt .* sarray) * transpose(Barray[2, :])


for tidx in 3:ntpoints - 1
	
	if mod(tidx,1000)==0
    	println("t-loop: \t $(round(Int, tidx/ntpoints * 100))%")
	end
 
    t = tarray[tidx]

    # Rk4 weights
    k1 = fvecLB(t, μarray[tidx, :])
    k2 = fvecLB(t + Δt/2, μarray[tidx, :] .+ Δt_finer/2 .* k1)
    k3 = fvecLB(t + Δt/2, μarray[tidx, :] .+ Δt_finer/2 .* k2)
    k4 = fvecLB(t + Δt, μarray[tidx, :] .+ Δt_finer .* k3)

    μarray[tidx+1, :] += μarray[tidx, :] .+ Δt_finer/6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4) + Δt_finer .* Barray[tidx, :]  # SEMIEXPLICIT
    #@assert minimum(μarray) > -1e-7 "minimum(μarray) = $(minimum(μarray)) at tidx = $tidx."
    # SEMIEXPLICIT >>>
    Bnext = zeros(ncells)
    #if tidx > 3
    if ncells == 1
        Bnext = A_finer \ (ℳ_finer*Barray[tidx, :] .+ 𝒩_finer*Barray[tidx-1, :] .+ sum(earray .* b4_finer) .* Barray[tidx-2, :] .+ Δt_finer .* Γ*μarray[tidx+1, :] .- Δt_finer .* [sum(earray[l] * exp(sarray[l]*Δt_finer) * He[l, j] for l in 1:2*M+1) for j in 1:ncells])
    else
        Bnext = A_finer \ (ℳ_finer*Barray[tidx, :] .+ 𝒩_finer*Barray[tidx-1, :] .+ sum(earray .* b4_finer) .* Barray[tidx-2, :] .+ Δt_finer .* Γ*μarray[tidx+1, :] .- Δt_finer .* [sum(earray[l] * exp(sarray[l]*Δt_finer) * He[l, j] for l in 1:2*M+1) for j in 1:ncells] .+ Δt_finer .* [sum(sum(wmatrix[j, k, l] * exp(sarray[l]*Δt_finer) * Hw[l, k] for l in 1:2*M+1) for k in filter!(x -> x≠j, collect(1:ncells))) for j in 1:ncells])
    end

    @assert maximum(abs.(imag.(Bnext))) < 1e-8 "maximum(abs.(imag.(Bnext))) = $(maximum(abs.(imag.(Bnext))))"
    Barray[tidx+1, :] .= real.(Bnext)
    # <<<

    if ncells > 1
        Hw .= Diagonal(exp.(Δt_finer .* sarray)) * Hw .+ b3_finer * transpose(Barray[tidx, :]) .+ b4_finer * transpose(Barray[tidx-1, :])
    end
    if tidx > 3
        He .= Diagonal(exp.(Δt_finer .* sarray)) * He .+ 1/Δt_finer .* b3_finer * transpose(Barray[tidx, :] .- Barray[tidx-1, :]) .+ 1/Δt_finer .* b4_finer * transpose(Barray[tidx-1, :] .- Barray[tidx-2, :])
    else  # tidx = 3
        He .= Diagonal(exp.(Δt_finer .* sarray)) * He .+ 1/Δt_finer .* b3_finer * transpose(Barray[3, :] .- Barray[2, :]) .+ b4_finer * transpose(dB2)
    end
end

writedlm("./datasets/v_array_freq$(freq)_$(n_on_fiber)line.txt", μarray)


