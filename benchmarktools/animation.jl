
function animation(sol; d::Real = 10, fps::Int = 30, truesol = nothing)

    # INPUTS
    #   - sol       : Solution of an HJ equation,
    #   - d         : Duration wanted for the animation,
    #   - fps       : Number of frame per secondes,
    #   - truesol   : if not nothing, it is the matrix of the true known solution

    @unpack Nt, T, space = sol.params  
    @unpack U = sol
    
    # set the min and maximum values of the plot, to scale the axis
    function scale_axis(V; _min = min(V...), _max = max(V...))
        minV = min(min(V...), _min)
        maxV = max(max(V...), _max)
        minV *= 1.05^sign(-minV)
        maxV *= 1.05^sign(maxV)
        [minV, maxV]
    end

    axis = scale_axis(U)

    # Total frame
    nbframe = fps * d

    # Space and time vector
    rescaltime = range(0,1,Nt+1)
    xLim = [space[begin], space[end]]

    anim = @animate for n in 1:nbframe

        t = n / (fps*d)

        # Computation of interpolations
        interpol = [linear_interpolation(rescaltime, U[i,:])(t) for i in eachindex(U[:,1])]

        # Create plots at rescaled time t
        pltsol = plot(space, interpol, label = sol.name, xlim = xLim, ylim= axis, legend=:bottomleft, lc=:black)

        # Plot of true sol
        if !isnothing(truesol)

        end
        pltsol
    end

    gif(anim, fps = fps)
end