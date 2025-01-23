module HamiltonJacobi1D

    using LinearAlgebra     # For matrix computations
    using UnPack            # Useful tool to unpack structure in Julia

    abstract type HJScheme end
     
    export HJParameters, HJProblem, HJSolution, HJEquation

    Base.@kwdef struct HJEquation
        H                 # Principal hamiltonian
        G  = nothing      # Second hamiltonian
        u0                # Initial Condition
        f  = nothing      # Source Term
    end

    include("parameters.jl")

    struct HJProblem{TS <: HJScheme}
        params::HJParameters        # Parameters
        equation::HJEquation        # Equation
        scheme::TS                  # Scheme
        name::String                # Name given to the problem
        function HJProblem(;T, Nt, L, Nx, H, u0, scheme, name = "")
            params = HJParameters(;T = T, Nt = Nt, L = L, Nx = Nx)
            equation = HJEquation(H = H, u0 = u0)
            new{typeof(scheme)}(params, equation, scheme, name)
        end
        function HJProblem(params::HJParameters, equation::HJEquation, scheme::HJScheme, name = "")
            new{typeof(scheme)}(params, equation, scheme, name)
        end
    end
    function HJProblem(prob::HJProblem; params = prob.params)
        HJProblem(params, prob.equation, prob.scheme, prob.name)
    end
    
    include("solution.jl")

    solve(prob::HJProblem) = solve(prob.params, prob.equation, prob.scheme, prob.name)
    
    function solve(params::HJParameters, equation::HJEquation, scheme::HJScheme, name::String = "")
        @unpack Nx, Nt, space = params
        @unpack u0 = equation
        # ALLOCATION FOR THE SOLUTION
        U  = zeros(Float64, Nx, Nt+1)
        # INITIALIZATION
        U[:,1]  .= u0.(space)
        # PERFORMSTEP (Currently it is all steps)
        performstep!(scheme, params, equation, U)
        # Make Solution
        HJSolution(params, equation, scheme, U, name)
    end

    export Upwind
    export solve
    include("schemes.jl")


end