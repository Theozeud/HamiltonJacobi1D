module HamiltonJacobi1D

    using LinearAlgebra     # For matrix computations
    using UnPack            # Useful tools to handle structure in Julia

    abstract type HJScheme end
     
    export HJParameters, HJProblem, HJSolution

    include("parameters.jl")

    struct HJProblem{TS <: HJScheme}
        params::HJParameters        # Parameters
        scheme::TS                  # Scheme
        name::String                # Name given to the problem
        function HJProblem(;T, Nt, L, Nx, H, u0, scheme, name = "")
            new{typeof(scheme)}(HJParameters(;T = T, Nt = Nt, L = L, Nx = Nx, H = H, u0 = u0), scheme, name)
        end
        function HJProblem(params::HJParameters, scheme::HJScheme, name = "")
            new{typeof(scheme)}(params, scheme, name)
        end
    end
    
    include("solution.jl")

    solve(params::HJParameters, scheme::HJScheme, name::String = "") = solve(Problem(params, scheme, name))
    solve(prob::HJProblem) = solve(prob.scheme, prob)

    export Upwind
    export solve
    include("schemes.jl")


end