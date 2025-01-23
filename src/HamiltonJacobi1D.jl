module HamiltonJacobi1D

    using LinearAlgebra     # For matrix computations
    using UnPack            # Useful tools to handle structure in Julia
     
    export Parameters, Problem

    include("parameters.jl")

    abstract type HJScheme end

    struct Problem{TS <: HJScheme}
        params::Parameters
        scheme::TS
        function Problem(;T, Nt, L, Nx, H, u0, scheme)
            new{typeof(scheme)}(Parameters(;T = T, Nt = Nt, L = L, Nx = Nx, H = H, u0 = u0), scheme)
        end
    end
    
    solve(scheme::HJScheme, ::Parameters) = @error "No solve function for the scheme $(typeof(scheme))"
    solve(prob::Problem) = solve(prob.scheme, prob.params)

    export Upwind
    export solve
    include("schemes.jl")


end