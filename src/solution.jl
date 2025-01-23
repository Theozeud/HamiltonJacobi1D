struct HJSolution
    params          # Parameters used
    equation        # Equation
    scheme          # Scheme used
    U               # Matrix of the solution
    name            # Name given to the simulation
    function HJSolution(params::HJParameters, equation::HJEquation, scheme::HJScheme, U, name)
        new(params, equation, scheme, U, name)
    end
end