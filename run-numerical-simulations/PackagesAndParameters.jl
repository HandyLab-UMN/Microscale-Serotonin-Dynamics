# Packages
using LinearAlgebra
using SpecialFunctions
using Distributions
using Dierckx
using Plots
using LaTeXStrings
using TickTock
using DelimitedFiles

include("rootfinder.jl")
include("keyFunctions.jl")

# parameters

# firing times
freq = 2      # standard options: 2, 16, 64
modelType = "nonlinear" # or "linear"

if (modelType != "nonlinear") & (modelType != "linear")
	error("Model type needs to be nonlinear or linear. You have: $(modelType)")
end

D = 1 # 3/30   #0.75          # diffusivity for u
σ = 0.04  #1.75 #1/5.8          # degradation rate for u

R0 = 0.3   # radius of compartment in μm
Len = 4.     # minimal separation of disks in μm, along fibers
Vol = 4/3*π*R0^3  # volume of ball with radius R0 in μm^3 = 10^(-15) L
εᵣ = R0/Len  #0.07  #0.05  # compartment radius
LF = 16/Len  # non-dim. parallel fiber separation
ν = -1/log(εᵣ)


# cell locations

# fibers with varicosities of radius εᵣ truly separated by x-distance 1
xc = []
yc = []
nfibers = 5
fiberindices = []
n_on_fiber = 1
xlength = LF*(nfibers - 1) + (nfibers - 1) * 2*εᵣ
ylength = (n_on_fiber - 1) + (n_on_fiber - 1) * 2*εᵣ
for vidx in 1:n_on_fiber
    global xc, yc
    xc = vcat(xc, [-xlength/2 + (j-1)/max(1, (nfibers-1)) * xlength for j in 1:nfibers])  # deterministic
    yc = vcat(yc, (-ylength/2 + (vidx-1)/max(1, (n_on_fiber-1)) * ylength)*ones(nfibers))
end
ncells = nfibers * n_on_fiber
cpos = [[xc[i], yc[i]] for i in 1:ncells]
# collecting x-directional varicosity indices of each fiber
for fiberidx in 1:nfibers
    append!(fiberindices, [[fiberidx + (j-1) * nfibers for j in 1:n_on_fiber]])
end
firingfibers = [1, 5]
nfiringfibers = length(firingfibers)
firingarray = zeros(Bool, nfiringfibers, ncells)  # marks all firing varicosity indices
firingvaric = zeros(Int, nfiringfibers * n_on_fiber)
for fireidx in 1:nfiringfibers
    varicidcs = [firingfibers[fireidx] + (j-1)*nfibers for j in 1:n_on_fiber]
    firingvaric[(fireidx-1)*n_on_fiber + 1:(fireidx-1)*n_on_fiber + n_on_fiber] .= varicidcs
    firingarray[fireidx, varicidcs] .= ones(Bool, n_on_fiber)
end


# universal constant
γe = Base.MathConstants.eulergamma

κ0vec = [1/(4*D) * exp(2*(γe)) for j in 1:ncells]               # for finite d₁: [1/(4*D) * exp(2*(γe - D/d₁vec[j])) for j in 1:ncells]


# time domain
kR = 300/Len^2 
tfinal = 1125.  # 60s ~ 1125
Δt = 0.03  # 0.03 # 0.01  # >> εᵣ^2 * κ0 * exp(-γₑ) = 0.0009 CHANGED
@assert 2*εᵣ^2*maximum(κ0vec)*exp(-γe) <  Δt &&  Δt < 1/2. "Δt = $Δt ≪ εᵣ^2*maximum(κ0vec)*exp(-γe) = $(εᵣ^2*maximum(κ0vec)*exp(-γe))"
#tarray = collect(0:Δt:tfinal) CHANGED
tfineness_factor = 5  # increasing the fineness of the time domain discretization
Δt_finer = Δt/tfineness_factor #/10
tarray = vcat([0, Δt], collect(2*Δt:Δt_finer:tfinal))
ntpoints = length(tarray)  # == tfinal/Δt + 1

freq_nondim = freq / kR 

Tfiber = []
nfirings = round(Int, (tfinal- 0.6*kR) * freq_nondim)  # no switching between frequencies

for fiberidx in firingfibers
    append!(Tfiber, [[0.6*kR + (k-1) * 1/freq_nondim for k in 1:nfirings]])  # no switching between frequencies
end

# Sel'kov reaction kinetics
# reaction kinetics: choosing Sel'kov intracellular kinetics
b = 250.  # release concentration kick in compartment
V = 10.67    # uptake rate from compartment
K = 7.5 * Vol # Michaelis Menten constant

std_kick = 0.03  # 0.015    # for smoothed δ(t - Tl)

fvec(t::Float64, μ::Array{Float64, 1}) = -V .* μ ./ (K .+ μ).*(modelType=="nonlinear") .+ -V .* μ ./K.*(modelType=="linear") .+ b/sqrt(2*π*std_kick^2) * sum(sum(exp(-(t - Tl)^2/(2*std_kick^2)) for Tl in Tfiber[fireidx]) .* firingarray[fireidx, :] for fireidx in 1:nfiringfibers)

fiberindex = firingfibers[1]
periodarray = Tfiber[fiberindex][2:end] .- Tfiber[fiberindex][1:end-1]
kick_averaged_array = [1/P for P in periodarray]
kick_averaged_array = vcat([0, 0, 0], kick_averaged_array, [0, 0, 0, 0])
tkicks = vcat([j/5 * Tfiber[fiberindex][1] for j in 1:3], Tfiber[fiberindex], [Tfiber[fiberindex][end] + j/5 * (tfinal - Tfiber[fiberindex][end]) for j in 3:5])
kick = Spline1D(tkicks, kick_averaged_array; k=1, s=0)  

ω = sqrt(σ/D)     # inverse of effective mobilization

# gauge function
ν = -1/log(εᵣ)

absxarray = zeros(ncells, ncells)
for i in 1:ncells
    for j in i+1:ncells
        # WARNING without "/ D", so not equal to aij
        absxarray[i, j] = sqrt((cpos[i][1] - cpos[j][1])^2 + (cpos[i][2] - cpos[j][2])^2)
        absxarray[j, i] = absxarray[i, j]
    end
end


# with different reaction rates for possibly each compartment
αvec = [2 * (1/ν + log(2*sqrt(D/σ)) - γe) for j in 1:ncells]    # for finite d₁: [2 * (1/ν + D/d₁vec[j] + log(2*sqrt(D/σ)) - γe) for j in 1:ncells]
C = Len^2 / (π*R0^2)                # = L² d₂/d₁
βvec = C * ones(ncells)     # for finite d₁: d₂vec ./ d₁vec
Γ = Diagonal(4*π*D .* βvec)


# discretization points of the Bromwich contour and weights to compute the heat
# kernel:
# G(t; a) = 1/(4πt) exp(-a²/(4t) - σt) = 1/(2πi) ∫ᵧ U(s) exp(st) ds
#Re(b3 w12) with U(s) = 1/(2π) K₀(a sqrt(s+σ))
δ = 1e-3  # time-interval start
T = 1e3   # time-interval end
m = 75    # number of discretization points: 2*m + 1
α = 0.8
β = 0.7
@assert 0 < α-β < α+β < π/2
θ = 0.95  # should be ∈ (0, 1)
a = acosh(2*T/(δ*(1-θ)*sin(α))) 
λ = 2*π*β*m*(1-θ) / (T*a)
h = a / m
sarray = [λ*(1 - sin(α + complex(0,1)*k*h)) for k in -m:m]
warray(absx::Real) = [λ*h/(2*π) * cos(α + complex(0,1)*k*h) * 2 * besselk(0, absx*sqrt((sarray[k+m+1]+σ)/D)) for k in -m:m]