struct HJSolution
    params          # Parameters used
    equation        # Equation
    scheme          # Scheme used
    U               # Matrix of the solution on the mesh
    UAna            # Matrix of the analytical solution on the mesh
    uAna            # Function of the analytical solution
    name            # Name given to the simulation
    function HJSolution(prob::HJProblem, U, name)
        @unpack time, space = prob.params
        @unpack u = prob.equation
        if !isnothing(u)
            UAna = zero(U)
            for i ∈ eachindex(time)
                for j ∈ eachindex(space)
                    UAna[j,i] = u(time[i], space[j])
                end
            end
            return new(prob.params, prob.equation, prob.scheme, U, UAna, u, name)
        else
            return new(prob.params, prob.equation, prob.scheme, U, nothing, nothing, name)
        end
    end
end