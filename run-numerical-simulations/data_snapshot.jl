# Find the snapshots; 
# Requires the other code to run for the arrays to exist but still takes a while since you have to step through time

include("PackagesAndParameters.jl")

tarray = readdlm("./datasets/tarray60s.txt")
μarray = readdlm("./datasets/mu_array_freq$(freq)_$(n_on_fiber)line_$(modelType).txt")
μarray = reshape(μarray, (length(tarray), ncells))
Barray = readdlm("./datasets/B_array_freq$(freq)_$(n_on_fiber)line_$(modelType).txt")
Barray = reshape(Barray, (length(tarray), ncells))

# steady-state # 614500
μeq = steadystate(freq_nondim, 4. .* sum(firingarray[fireidx, :] for fireidx in 1:nfiringfibers) .+ 1,modelType)  # [2., 1.5, 1., 1.5, 2.]
colors = [:magenta, :dodgerblue, :indianred, :orange, :purple, :orchid, :blue, :red, :darkorange1, :rebeccapurple]


# Plotting serotonin concentration on line connecting all varicosities
xarray = collect(-(5*εᵣ + 2.5*LF):0.05:5*εᵣ + 2.5*LF)  # so that not on cells
yarray = zeros(1)
# setting up an matrix with all weight arrays for all grid points and cells
wmatrix = zeros(Complex, length(xarray), length(yarray), ncells, 2*m+1)
for j in 1:ncells
    println("Computing exponential sum weights for the heat kernel: \t $(round(Int, j/ncells * 100))%")
    for xidx in eachindex(xarray)
        for yidx in eachindex(yarray)
            point = [xarray[xidx], yarray[yidx]]
            wmatrix[xidx, yidx, j, :] .= warray(norm(point .- cpos[j], 2))
        end
    end
end

Iarray = zeros(Complex, ncells, 2*m+1)  # 0 at t = 0
        
snapshot11 = true
snapshot21 = true
snapshot31 = true

snapshot12 = true
snapshot22 = true
snapshot32 = true

saveState = false

for tidx in 1:ntpoints - 1
    global Iarray, snapshot11, snapshot21, snapshot31, snapshot12, snapshot22, snapshot32, saveState
	
	if mod(tidx,1000)==0
    	println("line plot t-loop: \t $(round(Int, tidx/ntpoints * 100))%")
		println("$(tarray[tidx]) vs. $(Tfiber[fiberindex][2])")
		println("$(tarray[tidx]) vs. $(Tfiber[fiberindex][end-1])")
	end
	
    if tarray[tidx] < Tfiber[fiberindex][2]
        Uxy = zeros(length(xarray), length(yarray))
        for xidx in eachindex(xarray)
            for yidx in eachindex(yarray)
                incell = false
                for j in 1:ncells
                    if norm([xarray[xidx], yarray[yidx]] .- cpos[j], 2) < εᵣ
                        incell = true
                        Uxy[xidx, yidx] = μarray[tidx, j] * C
                    end
                end
                if incell == false
                    Uxyval = -1/(4*π*D) * sum(sum(wmatrix[xidx, yidx, j, l] * Iarray[j, l] for l in 1:2*m+1) for j in 1:ncells)
                    @assert abs(imag(Uxyval)) < 1e-8
                    Uxy[xidx, yidx] = real(Uxyval)
                    for idx in 1:nfiringfibers*n_on_fiber
                        dist_from_firing = norm([xarray[xidx], yarray[yidx]] .- cpos[firingvaric[idx]], 2)
                        if dist_from_firing < 1.9*εᵣ  # CHANGE DEPENDING ON FREQ, usually 10, 1.7
                            # choosing either global or local solution
                            Uxy[xidx, yidx] = min(Uxy[xidx, yidx], μarray[tidx, firingvaric[idx]] * C)  # Barray[tidx, firingvaric[idx]]/(2*π*D) * log(dist_from_firing/εᵣ) + C * μarray[tidx, firingvaric[idx]])  # μarray[tidx, firingvaric[idx]] * C
                        end
                    end
                end
            end
        end

        # snapshot plots 
        if Tfiber[fiberindex][2] > tarray[tidx] >= Tfiber[fiberindex][1]
            firingeventidx = 1  # length(Tfiber[fiberindex]) - 1
            # at firing 
            if snapshot11 == true
				stringName = "start1"
				saveState = true			
                snapshot11 = false
			elseif Tfiber[fiberindex][1] + 5/10 * 1/freq_nondim < tarray[tidx]  && snapshot21 == true
				stringName = "mid1"
				saveState = true
                snapshot21 = false
			elseif tarray[tidx] > Tfiber[fiberindex][1] + 8/10 * 1/freq_nondim && snapshot31 == true  # 8/10
				stringName = "end1"
				saveState = true
                snapshot31 = false
            end
        end
    end



	if tarray[tidx] >= Tfiber[fiberindex][end-1]  # && tarray[tidx] < Tfiber[fiberindex][2]  #&& tarray[tidx] > 19/20 * tfinal
        Uxy = zeros(length(xarray), length(yarray))
        for xidx in eachindex(xarray)
            for yidx in eachindex(yarray)
                incell = false
                for j in 1:ncells
                    if norm([xarray[xidx], yarray[yidx]] .- cpos[j], 2) < εᵣ
                        incell = true
                        Uxy[xidx, yidx] = μarray[tidx, j] * C
                    end
                end
                if incell == false
                    Uxyval = -1/(4*π*D) * sum(sum(wmatrix[xidx, yidx, j, l] * Iarray[j, l] for l in 1:2*m+1) for j in 1:ncells)
                    @assert abs(imag(Uxyval)) < 1e-8
                    Uxy[xidx, yidx] = real(Uxyval)
                    for idx in 1:nfiringfibers*n_on_fiber
                        dist_from_firing = norm([xarray[xidx], yarray[yidx]] .- cpos[firingvaric[idx]], 2)
                        if dist_from_firing < 1.9*εᵣ  # CHANGE DEPENDING ON FREQ, usually 10, 1.7
                            # choosing either global or local solution
                            Uxy[xidx, yidx] = min(Uxy[xidx, yidx], μarray[tidx, firingvaric[idx]] * C)  # Barray[tidx, firingvaric[idx]]/(2*π*D) * log(dist_from_firing/εᵣ) + C * μarray[tidx, firingvaric[idx]])  # μarray[tidx, firingvaric[idx]] * C
                        end
                    end
                end
            end
        end

        # snapshot plots 
        if tarray[tidx] >= Tfiber[fiberindex][end-1]  #+ kR / (2*64)  #REPETITION
            firingeventidx = 1  # length(Tfiber[fiberindex]) - 1
            # at firing 
            if snapshot12 == true
				stringName = "start2"
				saveState = true			
                snapshot12 = false
			elseif Tfiber[fiberindex][end-1] + 5/(10*freq_nondim) < tarray[tidx]  && snapshot22 == true  # + kR / (2*64)
				stringName = "mid2"
				saveState = true
                snapshot22 = false
            elseif tarray[tidx] > Tfiber[fiberindex][end-1] + 8/(10*freq_nondim) && snapshot32 == true # + kR / (2*64)
				stringName = "end2"
				saveState = true
                snapshot32 = false
            end
        end
    end
	
	if saveState
		writedlm("./datasets/Uxy$(freq)_$(n_on_fiber)line_$(modelType)_$(stringName).txt", Uxy[:, 1] ./ Len^2)
		writedlm("./datasets/scatterUEq$(freq)_$(n_on_fiber)line_$(modelType)_$(stringName).txt",  μeq .* C / Len^2)
		writedlm("./datasets/scatterU$(freq)_$(n_on_fiber)line_$(modelType)_$(stringName).txt",  μarray[tidx, :] * C / Len^2)
		
		saveState = false;
	end
    
    # updating
    for j in 1:ncells
        Iarray[j, :] .= [Iarray[j, l] * exp(sarray[l]*(tarray[tidx+1]-tarray[tidx])) for l in 1:2*m+1] + [Barray[tidx, j] * (exp(sarray[l]*(tarray[tidx+1]-tarray[tidx]))-1)/sarray[l] + (Barray[tidx+1, j]-Barray[tidx, j]) * (exp(sarray[l]*(tarray[tidx+1]-tarray[tidx]))-1-sarray[l]*(tarray[tidx+1]-tarray[tidx]))/(sarray[l]^2*(tarray[tidx+1]-tarray[tidx])) for l in 1:2*m+1]
    end
end

writedlm("./datasets/xVec_snapshot.txt", xarray .* Len)
writedlm("./datasets/varLocX_snapshot.txt", xc .* Len)

