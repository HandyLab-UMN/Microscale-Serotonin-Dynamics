# This file provides the code for finding the root of a function f: ℝᵈ → ℝ and also of f: ℝᵈ → ℝᵈ
# based on gradient descent to find starting points for the newton method and various other algorithms

using Optimization, ForwardDiff, Zygote
#using OptimizationOptimisers  # Or different optimizers: https://docs.sciml.ai/Optimization/stable/examples/rosenbrock/
using OptimizationOptimJL
using LinearAlgebra
using NonlinearSolve, SteadyStateDiffEq
using NLsolve

function getroot(func, x0)
    # returns one root of func
    # x0 is a less well-informed starting point guess

    f = OptimizationFunction(func, Optimization.AutoForwardDiff())
    #l1 = func(x0)
    prob = OptimizationProblem(f, x0)

    #optf = OptimizationFunction((u,_)->func(u), Optimization.AutoForwardDiff())
    #prob = OptimizationProblem(optf, x0)
    #sol = solve(prob, Adam(0.05), maxiters = 100, progress = false)

    #optf = OptimizationFunction((u,_)->func(u), Optimization.AutoZygote())
    #prob = OptimizationProblem(optf, x0)
    #sol = solve(prob, Adam(0.05), maxiters = 10, progress = true)

    ##optf = OptimizationFunction((u,_)->func(u), Optimization.AutoForwardDiff())
    #prob = OptimizationProblem(optf, x0)
    #sol = solve(prob, GradientDescent(), maxiters = 5, progress = false)

    # a few more Newton steps
    ##prob = OptimizationProblem(optf, x0)
    sol = solve(prob, Newton(), maxiters = 200)

    ##prob = OptimizationProblem(optf, sol.u)
    ##sol = solve(prob, BFGS(), maxiters = 200, progress = true)

    return sol.u
end

function getrootl1(funcl2, funcl1, x0)
    # returns one root of func
    # funcl2 is the l2 norm of the multidimensional function func and funcl1 is
    # the l1 norm
    # x0 is a less well-informed starting point guess

    # start with the more convex funcl2
    optf = OptimizationFunction((u,_)->funcl2(u), Optimization.AutoForwardDiff())
    prob = OptimizationProblem(optf, x0)
    sol = solve(prob, Newton(), maxiters = 200)
    # proceed with funcl1 for more precision
    optf = OptimizationFunction((u,_)->funcl1(u), Optimization.AutoForwardDiff())
    prob = OptimizationProblem(optf, sol.u)
    sol = solve(prob, Newton(), maxiters = 5)

    return sol.u
end

function getrootl1(funcvec, x0)
    # returns one root of funcvec: ℝᵈ → ℝᵈ
    # x0 is a less well-informed starting point guess

    # start with the more convex funcl2
    optf = OptimizationFunction((u,_)->norm(funcvec(u), 2), Optimization.AutoForwardDiff())
    prob = OptimizationProblem(optf, x0)
    sol = Optimization.solve(prob, Newton(), maxiters = 200)
    # proceed with funcl1 for more precision
    optf = OptimizationFunction((u,_)->norm(funcvec(u), 1), Optimization.AutoForwardDiff())
    prob = OptimizationProblem(optf, sol.u)
    sol = Optimization.solve(prob, Newton(), maxiters = 5)

    return sol.u
end

function getrootnograd(func, x0)
    #prob = OptimizationProblem((u,_)->func(u), x0)
    #sol = solve(prob, NelderMead())

    optf = OptimizationFunction((u,_)->func(u))
    prob = OptimizationProblem(optf, x0)
    sol = solve(prob, ParticleSwarm(n_particles = 100))

    return sol.u
end

# own implementation

function newton2D(rootfunc, startpnt::Tuple{Float64, Float64}, tol::Float64, stepcdiff2D::Tuple{Float64, Float64}) :: Tuple{Float64, Float64}
    # returns the root location of rootfunc(var1; var2) using Newton's method
    # and the centre difference for the Jacobian
    # When tol > |(var1ₙ, var2ₙ) - (var1ₙ₋₁, var2ₙ₋₁)|₂, (var1ᵢₙ, var2ₙ) is
    # returned
    # WARNING: this function only returns one root to each startpnt

    if norm(rootfunc(startpnt[1], startpnt[2]), 2) < 1e-8
        # a root is already found
        return startpnt
    end

    jacob = zeros(2,2)
    (Δvar1, Δvar2) = stepcdiff2D
    (var1bef, var2bef) = (Inf, Inf)
    (var1n, var2n) = startpnt
    idx = 0
    while norm([var1n-var1bef, var2n-var2bef]) >= tol
        idx = idx + 1
        println("Newton iteration ", idx)
        jacob[1,1] = (rootfunc(var1n+Δvar1, var2n)[1] - rootfunc(var1n-Δvar1, var2n)[1])/(2*Δvar1)
        jacob[1,2] = (rootfunc(var1n, var2n+Δvar2)[1] - rootfunc(var1n, var2n-Δvar2)[1])/(2*Δvar2)
        jacob[2,1] = (rootfunc(var1n+Δvar1, var2n)[2] - rootfunc(var1n-Δvar1, var2n)[2])/(2*Δvar1)
        jacob[2,2] = (rootfunc(var1n, var2n+Δvar2)[2] - rootfunc(var1n, var2n-Δvar2)[2])/(2*Δvar2)
        # check invertibility
        if abs(det(jacob)) < 1e-14
            error("D(rootfunc) is not invertible for (var1, var2)=($var1n, $var2n)")
        end

        nextpnt = [var1n; var2n] - jacob \ rootfunc(var1n, var2n)

       # update approximation
       var1bef = var1n
       var2bef = var2n
       var1n = nextpnt[1]
       var2n = nextpnt[2]
       println("rootfunc(var1n, var2n) = ", rootfunc(var1n, var2n))
    end

    return var1n, var2n
end

function getnewtstartpnts(rootfunc, var1arraynewt::Array{Float64, 1}, var2arraynewt::Array{Float64, 1}, thresh::Float64) :: Array{Tuple{Float64,Float64}, 1}
    # finds starting points near roots for Newton's method
    # WARNING: Take var1arraynewt and var2arraynewt fine enough so that roots
    # are separated by at least two grid points

    rootfuncmatrix = zeros(length(var1arraynewt), length(var2arraynewt))
    println("Finding start points... ")
    for (var1idx, var1val) in enumerate(var1arraynewt)
        println("\t start point: ", round(Int,var1idx/length(var1arraynewt)*100), "\t %")
        for (var2idx, var2val) in enumerate(var2arraynewt)
            fval = rootfunc(var1val, var2val)
            rootfuncmatrix[var1idx, var2idx] = fval[1]^2 + fval[2]^2
        end
    end
    #Hopfrootplot = Plots.heatmap(var1arraynewt, var2arraynewt, rootfuncmatrix', title="fval[1]^2 + fval[2]^2", xlabel="var1", ylabel="var2", clims=(-thresh, thresh), color=:seismic, dpi=300)
    #png(Hopfrootplot, "HopfrootplotBravaisRM.png")

    indcs = findall(x -> x<thresh, rootfuncmatrix)
    @assert length(indcs) > 1

    # extract the starting points to each found root cluster/box
    notyetcaptured = true
    i1start = indcs[1][1]
    i2start = indcs[1][2]
    i1end = 0.
    i2end = 0.
    startpnts = []
    for i in 1:length(indcs)-1
        if indcs[i+1][1] - indcs[i][1] != 1 && notyetcaptured
            # jumped in first index either to next root box or next second index
            i1end = indcs[i][1]
            var1start = var1arraynewt[i1start]
            var1end = var1arraynewt[i1end]
            var2start = var2arraynewt[i2start]
            push!(startpnts, (mid(var1start..var1end), var2start))  # only taking start point in second variable not the centre of the box
            # update root box parameters
            i1start = indcs[i+1][1]
            if indcs[i+1][2] - indcs[i][2] == 1
                # jumped to next second index in same root boxes
                notyetcaptured = false
            end
        end
        if indcs[i+1][2] - indcs[i][2] > 1
            # arrived at next uncaptured root box having moved along second indices
            # update root box parameters
            notyetcaptured = true
            i1start = indcs[i+1][1]
            i2start = indcs[i+1][2]
        end
    end
    println(startpnts)

    return startpnts
end

function newton_nD(rootfunc, startpnt::Array{Float64, 1}, tol::Float64, stepcdiff_nD::Array{Float64, 1}) :: Array{Float64, 1}
    # WARNING Not yet written
    # returns the root location of rootfunc(var1, var2, ..., varn) using Newton's method
    # and the centre difference for the Jacobian
    # When tol > |(var1ₙ, ..., varnₙ) - (var1ₙ₋₁, ..., varnₙ₋₁)|₂, (var1ᵢₙ, var2ₙ) is
    # returned
    # WARNING: this function only returns one root to the startpnt

    if norm(rootfunc(startpnt...), 2) < 1e-8
        # a root is already found
        return startpnt
    end

    dim = length(startpnt)

    jacob = zeros(2, dim)
    Δpnt = stepcdiff_nD
    pntbef = Inf * ones(dim)
    pntnow = startpnt
    idx = 0
    while norm(pntnow-pntbef, 2) >= tol
        idx = idx + 1
        #println("Newton iteration ", idx)
        for j in 1:dim
            ej = zeros(dim)
            ej[j] = 1
            jacob[:, j] = (rootfunc(pntnow + stepcdiff_nD[j] .* ej ...) .- rootfunc(pntnow - stepcdiff_nD[j] .* ej ...)) ./ (2*stepcdiff_nD[j])
        end

        # check invertibility
        uJ, sJ, vJ = svd(jacob)
        if minimum(abs.(sJ)) < 1e-14  # one of the singular values has |...| = 0
            error("D(rootfunc) is not (pseudo)-invertible for point = $pntnow.")
        end
        
        # Newton-Raphson step
        nextpnt = pntnow - jacob \ rootfunc(pntnow...)  # pseudo-inverse pinv(jacob)

       # update approximation
       pntbef = pntnow
       pntnow = abs.(nextpnt)  # brute-forcely not allowing negative values
       #println("rootfunc($pntnow) = ", rootfunc(pntnow...))
    end

    return pntnow
end

# other well-performing algorithms

function getroot_steadystate(rootfunc, u0::Array{Float64, 1}) ::Array{Float64, 1}
    f(u, p, t) = rootfunc(u)
    prob = SteadyStateProblem(f, u0)
    sol = solve(prob, SSRootfind())

    return sol.u
end

function getroot_NL(rootfunc, u0)
    sol = NLsolve.nlsolve(rootfunc, u0)

    return sol.zero
end


#rosenbrock(x) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
#x0 = zeros(4)
#_p = [1.0, 100.0]

#println("")
#startpnt = getrootl1(rosenbrock, x0)
#println("rosenbrock = ", rosenbrock(startpnt))
