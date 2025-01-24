using HamiltonJacobi1D

include("../benchmarktools/setup.jl")


equation = minus_bump_equation(; domain = (-1,1))

problem = HJProblem(;
                    T           = 0.125,                            # Final Time
                    Nt          = 150,                         # Number of timesteps
                    Nx          = 500,                          # Nx space points
                    equation    = equation,                     # Equation to solve
                    scheme      = SemiLagrangian(),                  # Scheme
                    name        = "SemiLagrangian")                      # Abritrary Name of the problem

@time sol = solve(problem)