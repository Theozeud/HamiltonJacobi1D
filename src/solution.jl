struct HJSolution
    params          # Parameters used
    scheme          # Scheme used
    U               # Matrix of the solution
    name            # Name given to the simulation
    function HJSolution(prob::HJProblem, U)
        new(prob.params, prob.scheme, U, prob.name)
    end
end