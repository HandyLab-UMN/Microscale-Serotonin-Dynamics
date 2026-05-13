# Finds the lower bound of the steady state for the full ECS space (i.e., 2D domain)

include("PackagesAndParameters.jl")

function Bsteadystate(μeq::Array{Float64, 1}, ncells::Int) :: Array{Float64, 1}
    # returns the presumable steady-state close to guess μinit of the dynamics for time-averaged firing with a single frequency freq_nondim
    Me = zeros(ncells, ncells)
    Me .= Diagonal(-αvec)
    for j in 1:ncells-1
        for k in j+1:ncells
            Me[j, k] = -2*besselk(0, sqrt(σ/D)*absxarray[j, k])
            Me[k, j] = Me[j, k]
        end
    end
    Beq = Me \ (Γ * μeq)

    return Beq  
end

μeq = steadystate(freq_nondim, 4. .* sum(firingarray[fireidx, :] for fireidx in 1:nfiringfibers) .+ 1, modelType)  # [2., 1.5, 1., 1.5, 2.]
Beq = Bsteadystate(μeq, ncells)

# Plotting lower-bound time-averaged extracellular steady-state Ue(t, x)
xarray = collect(-15.:0.02:15.) 
yarray = collect(-15.:0.02:15.)
nxpnts = length(xarray)
nypnts = length(yarray)
Uarray = zeros(nxpnts, nypnts)
for xidx in 1:nxpnts
    for yidx in 1:nypnts
        in_nbhd = false
        for k in 1:ncells
            if norm([xarray[xidx], yarray[yidx]] .- cpos[k], 2) < εᵣ
                Uarray[xidx, yidx] = μeq[k] .* C 
                in_nbhd = true
            end
        end
        if in_nbhd == false
            Uarray[xidx, yidx] = - sum(Beq[j]/(2*π*D) * besselk(0, sqrt(σ/D) * norm([xarray[xidx], yarray[yidx]] .- cpos[j])) for j in 1:ncells)
        end
    end
end

writedlm("./datasets/xarray2D.txt", xarray.* Len)
writedlm("./datasets/yarray2D.txt", yarray.* Len)
writedlm("./datasets/U2D_freq$(freq)_$(n_on_fiber)lines_$(modelType).txt", Uarray' ./ Len^2)

# difference to single line (lower line for two (even) lines, middle line for three (odd) lines)
xc1 = zeros(nfibers)
yc1 = zeros(nfibers)
if mod(n_on_fiber, 2) == 0
    # take lower line
    xc1 .= xc[1:nfibers]
    yc1 .= yc[1:nfibers]
	
	filenametosave = "./datasets/U2D_freq$(freq)_$(n_on_fiber)lines_$(modelType)_even.txt"
else 
    # take middle line
    idxstart = Int(n_on_fiber/2 - 0.5)
    xc1 .= xc[idxstart*nfibers + 1:(idxstart+1)*nfibers]
    yc1 .= yc[idxstart*nfibers + 1:(idxstart+1)*nfibers]
	
	filenametosave = "./datasets/U2D_freq$(freq)_$(n_on_fiber)lines_$(modelType)_odd.txt"
end

n_on_fiber = 1
ncells = nfibers * n_on_fiber
cpos = [[xc1[i], yc1[i]] for i in 1:ncells]

firingfibers = [1, 5]
nfiringfibers = length(firingfibers)
firingarray = zeros(Bool, nfiringfibers, ncells)  # marks all firing varicosity indices
for fireidx in 1:nfiringfibers
    varicidcs = [firingfibers[fireidx] + (j-1)*nfibers for j in 1:n_on_fiber]
    firingarray[fireidx, varicidcs] .= ones(Bool, n_on_fiber)
end

αvec = [2 * (1/ν + log(2*sqrt(D/σ)) - γe) for j in 1:ncells]    # for finite d₁: [2 * (1/ν + D/d₁vec[j] + log(2*sqrt(D/σ)) - γe) for j in 1:ncells]
βvec = C * ones(ncells)     # for finite d₁: d₂vec ./ d₁vec
Γ = Diagonal(4*π*D .* βvec)


# Plotting lower-bound time-averaged extracellular steady-state Ue(t, x)
μeq = steadystate(freq_nondim, 4. .* sum(firingarray[fireidx, :] for fireidx in 1:nfiringfibers) .+ 1,modelType)  # [2., 1.5, 1., 1.5, 2.]
Beq_single = Bsteadystate(μeq, ncells)

U_single_array = zeros(nxpnts, nypnts)
for xidx in 1:nxpnts
    for yidx in 1:nypnts
        in_nbhd = false
        for k in 1:ncells
            if norm([xarray[xidx], yarray[yidx]] .- cpos[k], 2) < εᵣ
                U_single_array[xidx, yidx] = μeq[k] .* C 
                in_nbhd = true
            end
        end
        if in_nbhd == false
            U_single_array[xidx, yidx] = - sum(Beq_single[j]/(2*π*D) * besselk(0, sqrt(σ/D) * norm([xarray[xidx], yarray[yidx]] .- cpos[j])) for j in 1:ncells)
        end
    end
end

writedlm(filenametosave, U_single_array'./ Len^2)


