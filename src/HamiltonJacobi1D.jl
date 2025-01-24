module HamiltonJacobi1D

    using LinearAlgebra     # For matrix computations
    using UnPack            # Useful tool to unpack structure in Julia

    abstract type HJScheme end
     
    export HJParameters, HJProblem, HJSolution, HJEquation

    Base.@kwdef struct HJEquation
        domain::Tuple{Real,Real}    # Spatial domain of the equation
        H                           # Principal hamiltonian
        G  = nothing                # Second hamiltonian
        f  = nothing                # Source Term
        u0                          # Initial Condition
        u  = nothing                # Analytical solution if it is known
    end

    include("parameters.jl")

    struct HJProblem{TS <: HJScheme}
        params::HJParameters        # Parameters
        equation::HJEquation        # Equation
        scheme::TS                  # Scheme
        name::String                # Name given to the problem
        function HJProblem(params::HJParameters, equation::HJEquation, scheme::HJScheme, name = "")
            new{typeof(scheme)}(params, equation, scheme, name)
        end
        function HJProblem(;T, Nt, Nx, equation::HJEquation, scheme::HJScheme, name = "")
            params = HJParameters(;T = T, Nt = Nt, Linf = equation.domain[1], Lsup = equation.domain[2], Nx = Nx)
            new{typeof(scheme)}(params, equation, scheme, name)
        end
    end
    function HJProblem(prob::HJProblem; params::HJParameters = prob.params, 
        scheme::HJScheme = prob.scheme, name::String = prob.name)
        HJProblem(params, prob.equation,scheme, name)
    end
    
    include("solution.jl")

    solve(params::HJParameters, equation::HJEquation, scheme::HJScheme, name::String = "") =
        solve(HJProblem(params, equation, scheme, name))
    
    function solve(prob::HJProblem)
        @unpack params, equation, scheme, name = prob
        @unpack Nx, Nt, space = params
        @unpack u0 = equation
        # ALLOCATION FOR THE SOLUTION
        U  = zeros(Float64, Nx, Nt+1)
        # INITIALIZATION
        U[:,1]  .= u0.(space)
        # PERFORMSTEP (Currently it is all steps)
        performstep!(scheme, params, equation, U)
        # Make Solution
        HJSolution(prob, U, name)
    end

    export Upwind
    export solve
    include("schemes.jl")


end