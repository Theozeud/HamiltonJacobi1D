using HamiltonJacobi1D

include("../benchmarktools/setup.jl")

problem = HJProblem(;
                    T           = 1,                            # Final Time
                    Nt          = 10000,                        # Number of timesteps
                    L           = 1,                            # Domain is [0,L]
                    Nx          = 500,                          # Nx space points
                    H           = x -> abs(x)^2/2,              # Hamiltonian
                    u0          = x -> max(1-16*(x-0.5)^2,0),   # Initial condition
                    scheme      = Upwind(0.0))                  # Scheme
                    name        = "Upwind"                      # Abritrary Name of the problem

@time sol = solve(problem)