struct HJParameters
    T               # Final Time
    Nt              # Number of time iterations
    Δt              # Timestep
    time            # Timespan
    Linf            # Space domain is [Linf,Lsup]
    Lsup            # 
    Nx              # Size of the grids in real space
    Δx              # Spacestep
    space           # Space domain
    function HJParameters(;T, Nt, Linf, Lsup, Nx)
        Δt      = T/Nt
        time    = range(0,T,Nt+1)
        Δx      = (Lsup - Linf)/(Nx-1)
        space   = range(Linf,Lsup,Nx)
        new(T, Nt, Δt, time, Linf, Lsup, Nx, Δx, space)
    end
    #=
    function HJParameters(;T, Δt, Linf, Lsup, Δx)
        Nt      = Int(floor(T/Δt))
        _T      = Nt * Δt
        time    = range(0,T,Nt+1)
        Nx      = Int(floor((Lsup - Linf)/Δx))+1
        _Lsup   = (Nx-1) * Δx + Linf
        space   = range(Linf,Lsup,Nx)
        new(_T, Nt, Δt, time, Linf, _Lsup, Nx, Δx, space)
    end
    =#
end

function HJParameters(params; T = params.T, Nt = params.Nt, Linf = params.Linf, Lsup = params.Lsup, Nx = params.Nx)
    HJParameters(;T = T, Nt = Nt, Linf = Linf, Lsup = Lsup, Nx = Nx)
end

#=
function HJParameters(params; T = params.T, Δt = params.Δt, Linf = params.Linf, Lsup = params.Lsup, Δx = params.Δx)
    HJParameters(;T = T, Δt = Δt, Linf = Linf, Lsup = Lsup, Δx = Δx)
end
=#
function Base.show(io::IO, params::HJParameters)
    printstyled(io, "Parameters :"; bold = true)
    println(io,"")
    println(io, "Time = [0,$(params.T)]")
    println(io, "Nt = $(params.Nt)")
    println(io, "Δt = $(params.Δt)")
    println(io, "Space = [$(params.space[begin]),$(params.space[end])]")
    println(io, "Nx = $(params.Nx)")
    println(io, "Δx = $(params.Δx)")
end