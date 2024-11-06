"""

we're solving the following problem

∂ₜ u + div f(u,p,t) = s(u,∇u,p,t)

f(u,t,i,j) where i is the line of f and j is the collumn of f

"""

function f(u,p,t,i,j)
    h,qx,qy = u
    g = 9.81
    if i==1
        if j==1
            return qx
        else
            return qy
        end
    elseif i==2
        if j==1
            return qx^2/h+g*h^2/2
        else
            return qx*qy/h
        end
    else
        if j==1
            return qx*qy/h
        else
            return qy^2/h + g*h^2/2
        end
    end
end

function s(u,dux,duy,p,t,i)
    if i==1
        return 0.0
    elseif i==2
        return 0.0
    else
        return 0.0
    end
end

# we also need to define the eigen values of f called λ
@fastmath function λ(u,p,t,i,j)
    h,qx,qy = u
    g = 9.81
    if i==1
        if j==1
            return qx/h
        else
            return qy/h
        end
    elseif i==2
        if j==1
            return qx/h+sqrt(g*h)
        else
            return qy/h+sqrt(g*h)
        end
    else
        if j==1
            return qx/h-sqrt(g*h)
        else
            return qy/h-sqrt(g*h)
        end
    end
end

# About the initial condition now

function u0(p,i)
    if i==1
        return 1.0
    elseif i==2
        return 0.0
    else
        return 0.0
    end
end

# Finally the boundary conditions, we will use Dirichlet conditions here, however, take care of hyperbolic restrictions when using this.
# Also, for neumann and periodic conditions, we will get a u inside those function, for neumann the value that should be entered is the one at the symetric position in the element, for periodic, use the value at the other side of the mesh
# Final note, the point entering here should be a Point{} in the real domain, so that it can be used in general mesh

function uleft(u,p,t,i)
    if i==1
        return u[1] # neumann
    elseif i==2
        return u[i]#sinpi(4t)^2*exp(-50*(p.y-1.0)^2) # enter a wave function in the x direction
    else
        return u[i]
    end
end

function uright(u,p,t,i)
    return u[i] # let's neumann everywhere here
end

function ubot(u,p,t,i)
    if i==1
        return u[1] # neumann
    elseif i==2
        return u[i]#sinpi(4t)^2*exp(-50*(p.y-1.0)^2) # enter a wave function in the x direction
    else
        return   sinpi(4t)^2*exp(-50*(p.x-1.0)^2)#
    end
end

function utop(u,p,t,i)
    return u[i] # let's neumann everywhere here
end

function callback(df,cnt)
    return false
end