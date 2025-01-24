# SET OF EQUATION WITH KNOWN ANALYTICAL SOLUTION

# Lipschitz initial condition with bounded support
# Quadratic hamiltonian
function bump_equation(;domain::Tuple{<:Real,<:Real})
    Linf   = domain[1]
    Lsup   = domain[2]
    Lmid   = (Linf + Lsup)/2 
    ΔL     = Lsup - Linf
    H(x)   = abs(x)
    u0(x)  = max(ΔL^2-16*(x-Lmid)^2,0)
    function u(t,x)
        if t < 1/32
            if Linf + ΔL/4 ≤ x ≤ Lmid + ΔL/4
                return 16*(x - Lmid)^2/(32*t - 1) + ΔL^2
            else
                return zero(x)
            end
        else
            if Linf + ΔL/4 ≤ x ≤ Lmid + ΔL/4
                return 16*(abs(x-Lmid) - ΔL/4)^2/(32*t - 1)
            else
                return zero(x)
            end
        end
    end
    return HJEquation(domain = domain, H = H, u0 = u0, u = u)
end

# Lipschitz and semi-concave initial condition with bounded support 
# Quadratic hamiltonian
function minus_bump_equation(;domain::Tuple{<:Real,<:Real})
    Linf   = domain[1]
    Lsup   = domain[2]
    Lmid   = (Linf + Lsup)/2 
    ΔL     = Lsup - Linf
    H(x)   = abs(x)/2
    u0(x)  = -max(ΔL^2-16*(x-Lmid)^2,0)
    u(t,x) = min(abs(x-Lmid)^2/(t+1/16)-ΔL^2,0)
    return HJEquation(domain = domain, H = H, u0 = u0, u = u)
end



