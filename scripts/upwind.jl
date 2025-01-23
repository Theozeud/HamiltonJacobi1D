using HamiltonJacobi1D

include("benchmarktools.jl")

problem = Problem(;
                    T           = 2.5,                          # Final Time
                    Nt          = 100,                          # Number of timesteps
                    L           = 1,                            # Domain is [0,L]
                    Nx          = 1000,                          # Nx space points
                    H           = x -> abs(x)^2,              # Hamiltonian
                    u0          = x -> max(1-16*(x-0.5)^2,0),   # Initial condition
                    scheme      = Upwind(0.0))                  # Scheme

@time solve(problem)