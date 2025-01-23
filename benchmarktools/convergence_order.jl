struct ConvergenceResult
    ΔxΔt    # Couple of time/space steps
    cfl     # CFL Number Δt / Δx keeping constant
    Error   # Vector of errors
end

function conv(problems::Vector{Problem}; Δx, cfl, truesol)
    Errors = Dict()
    ΔxΔt = [(Δx[i],cfl * Δx[i]) for i ∈ eachindex(Δx)]
    for problem ∈ problems
        @unpack name = problem
        println(name)
        Error = zeros(T, length(Δx))
        @inbounds for i ∈ eachindex(Δx)
            newparams = HJParameters(problem.params; Δx = ΔxΔt[i][1], Δt = ΔxΔt[i][2])
            newprob = Problem(problem; params = newparams)
            @time "(Δx,Δt) = $(ΔxΔt[i])" sol = solve(newprob)
            TRUE_SOL = zero(U)
            @unpack time, space = newparams
            for i ∈ axes(U,2)
                TRUE_SOL[:,i] = [truesol(time[i], space[i]) for i ∈ eachindex(space)]
            end
            Error[i] = max((sol.U .- TRUE_SOL)...)
        end
        Errors[name] = Error
    end
    ConvergenceResult(ΔxΔt, cfl, Errors)
end