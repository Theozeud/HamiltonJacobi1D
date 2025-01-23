# UPWIND SCHEME

struct Upwind <: HJScheme 
    α0    # Parameter 
end

function numericalHamiltonian(scheme::Upwind, H, x, y)
    @unpack α0 = scheme
    if x ≥ α0 && y ≥ α0 
        return H(x)
    elseif x ≥ α0 && y < α0
        return  H(x) +  H(y)
    elseif x < α0 && y ≥ α0
        return H(α0)
    else
        return H(y)
    end
end


function solve(scheme::Upwind, params::Parameters)
    # DISCRETIZATION PARAMETERS
    @unpack Nt, Δt, Nx, Δx, H, U0 = params
    # ALLOCATION
    U       = zeros(Float64, Nx, Nt+1)
    # STORE INITIAL STATE
    U[:,1]  .= U0
    # LOOP FOR TIME ITERATIONS
    for n in 1:Nt
        U[1,n+1] = U[1,n]  - Δt * numericalHamiltonian(scheme, H, (U[1,n] - U[end,n]) / Δx , (U[2,n] - U[1,n]) / Δx)
        for j ∈ 2:Nx-1
            U[j,n+1] = U[j,n]  - Δt * numericalHamiltonian(scheme, H, (U[j,n] - U[j-1,n]) / Δx , (U[j+1,n] - U[j,n]) / Δx)
        end
        U[Nx,n+1] = U[Nx,n]  - Δt * numericalHamiltonian(scheme, H, (U[Nx,n] - U[Nx-1,n]) / Δx , (U[1,n] - U[Nx,n]) / Δx)
    end
    U
end