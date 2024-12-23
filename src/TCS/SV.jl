"""

we're solving the following problem

∂ₜ u + div f(u,p,t) = s(u,∇u,p,t)

f(u,t,i,j) where i is the line of f and j is the collumn of f

"""

bat(p) = exp(-50*(p.x-1.0)^2-50*(p.y-1.0)^2)

batpx(p) = ForwardDiff.derivative(x->bat(Point(x,p.y)),p.x)
batpy(p) = ForwardDiff.derivative(y->bat(Point(p.x,y)),p.y)
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
    h,qx,qy = u
    g=9.81
    if i==1
        return 0.0
    elseif i==2
        return -g*h*batpx(p)
    else
        return -g*h*batpy(p)
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
        return 2.0-bat(p)
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
        return 10*sinpi(4t)^2*exp(-50*(p.x-1.0)^2)#
    end
end

function utop(u,p,t,i)
    return u[i] # let's neumann everywhere here
end

function callback(df,cnt)
    ndf = @chain df begin
        @mutate(Height = u1+bat(!!Point(x,y)))
        @mutate(height = u1,qx = u2,qy = u3, bati = !!bat(!!Point(x,y)))
        @mutate(ux = qx/height,uy = qy/height)
        @select(x,y,Height,height,ux,uy,qx,qy,bati,time = t)
    end
    CSV.write("./datasCyl/data$cnt.csv",ndf)
    return false
end