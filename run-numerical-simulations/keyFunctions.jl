
# function fvec_(t::Float64, μ::Array{Float64, 1}) :: Array{Float64, 1}
#     return -V .* μ ./ (K .+ μ) .+ b * kick(t) .* sum(firingarray[fireidx, :] for fireidx in 1:nfiringfibers)
# end

function get_spikeheight(μeq::Array{Float64, 1}, std_kick::Float64) :: Vector{Float64}
    # returns the spike height seen from E∞[u]

    Me = zeros(ncells, ncells)
    Me .= Diagonal(-αvec)
    for j in 1:ncells-1
        for k in j+1:ncells
            Me[j, k] = -2*besselk(0, sqrt(σ/D)*absxarray[j, k])
            Me[k, j] = Me[j, k]
        end
    end
    Beq = Me \ (Γ * μeq)
    
    ΔT = 4*std_kick
    IGheat(x1, x2) = expint(1, norm(x1 .- x2, 2)^2/(4*D * ΔT)) * (1 + norm(x1 .- x2, 2)^2/(4*D/σ)) - exp(-norm(x1 .- x2, 2)^2/(4*D * ΔT)) * σ * ΔT
    Matri = zeros(ncells, ncells)
    Matri .= Diagonal(expint(1, σ * ΔT) .- (exp(-σ * ΔT) - 1) / (σ * ΔT).- αvec)
    for j in 1:ncells
        for k in j+1:ncells
            Matri[j,k] = -1/2*IGheat([xc[j], yc[j]], [xc[k], yc[k]])
            Matri[k,j] = Matri[j,k]
        end
    end
    
    NormDist = Normal(0, 1)
    uvec = (I - 2*ΔT/3 * inv(Matri)*Γ) \ ((2*cdf(NormDist, 2)-1  + freq_nondim*ΔT) * [b, 0, 0, 0, b] .- V.*μeq ./ (K .+ μeq) * ΔT  .+ ΔT * Beq)
    return uvec
end

function steadystate(freq_nondim::Float64, μinit::Array{Float64, 1},modelType::String) :: Array{Float64, 1}
    # returns the presumable steady-state close to guess μinit of the dynamics for time-averaged firing with a single frequency freq_nondim
    
    kick_averaged1 = freq_nondim * 1  # not accurate: freq_nondim * (cdf(Normal(0, std_kick), 2/freq_nondim * 1/2) - cdf(Normal(0, std_kick), -2/freq_nondim * 1/2))
    # differently firing varicosities
    #kick_averaged2 = freq_nondim2 * 1  # not accurate: freq_nondim2 * (cdf(Normal(0, std_kick), 2/freq_nondim2 * 1/2) - cdf(Normal(0, std_kick), -2/freq_nondim2 * 1/2))

    fvecSS(μ) = -V .* μ ./ (K .+ μ).*(modelType=="nonlinear") .+ -V .* μ ./K.*(modelType=="linear") .+ b * kick_averaged1 .* sum(firingarray[fireidx, :] for fireidx in 1:nfiringfibers)

    Me = zeros(ncells, ncells)
    Me .= Diagonal(-αvec)
    for j in 1:ncells-1
        for k in j+1:ncells
            Me[j, k] = -2*besselk(0, sqrt(σ/D)*absxarray[j, k])
            Me[k, j] = Me[j, k]
        end
    end
    rootfunc(μe) = Me * fvecSS(μe) .+ Γ * μe  # do not use any type info like :: Array
    
    return getroot_NL(rootfunc, μinit)
end