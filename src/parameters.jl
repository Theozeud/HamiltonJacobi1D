struct HJParameters
    T               # Final Time
    Nt              # Number of time iterations
    Δt              # Timestep
    time            # Timespan
    L               # Space domain is [0,L]
    Nx              # Size of the grids in real space
    Δx              # Spacestep
    space           # Space domain
    function HJParameters(;T, Nt, L, Nx)
        Δt      = T/Nt
        time    = range(0,T,Nt+1)
        Δx      = L/(Nx-1)
        space   = range(0,L,Nx)
        new(T, Nt, Δt, time, L, Nx, Δx, space)
    end
end


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