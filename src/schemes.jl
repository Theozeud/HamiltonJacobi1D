# UPWIND SCHEME

struct Upwind <: HJScheme 
    α0    # Parameter 
end

function numericalHamiltonian(scheme::Upwind, H, x, y)
    @unpack α0 = scheme
    if x ≥ α0 && y ≥ α0 
        return H(x)
    elseif x ≥ α0 && y < α0
        return  H(x) +  H(y) - H(α0)
    elseif x < α0 && y ≥ α0
        return H(α0)
    else
        return H(y)
    end
end


function performstep!(scheme::Upwind, params::HJParameters, equation::HJEquation, U)
    # DISCRETIZATION PARAMETERS
    @unpack Nt, Δt, Nx, Δx = params
    @unpack H = equation
    # LOOP FOR TIME ITERATIONS
    for n in 1:Nt
        U[1,n+1] = U[1,n]  - Δt * numericalHamiltonian(scheme, H, (U[1,n] - U[end,n]) / Δx , (U[2,n] - U[1,n]) / Δx)
        for j ∈ 2:Nx-1
            U[j,n+1] = U[j,n]  - Δt * numericalHamiltonian(scheme, H, (U[j,n] - U[j-1,n]) / Δx , (U[j+1,n] - U[j,n]) / Δx)
        end
        U[Nx,n+1] = U[Nx,n]  - Δt * numericalHamiltonian(scheme, H, (U[Nx,n] - U[Nx-1,n]) / Δx , (U[1,n] - U[Nx,n]) / Δx)
    end
    nothing
end

# SEMI-LAGRANGIAN

struct SemiLagrangian <: HJScheme end

function performstep!(::SemiLagrangian, params::HJParameters, equation::HJEquation, U)
    # DISCRETIZATION PARAMETERS
    @unpack Nt, Δt, Nx, Δx = params
    @unpack H = equation
    # LOOP FOR TIME ITERATIONS
    for n in 1:Nt
        U[1,n+1] = U[1,n]  - Δt * max( H(0), H((U[2,n] - U[1,n])/Δx), H(-(U[1,n] - U[Nx,n])/Δx))
        for j ∈ 2:Nx-1
            U[j,n+1] = U[j,n]  - Δt * max( H(0), H((U[j+1,n] - U[j,n])/Δx), H(-(U[j,n] - U[j-1,n])/Δx))
        end
        U[Nx,n+1] = U[Nx,n]  - Δt * max( H(0), H((U[1,n] - U[Nx,n])/Δx), H(-(U[Nx,n] - U[Nx-1,n])/Δx))
    end
    nothing
end
