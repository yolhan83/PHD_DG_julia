"""

we're solving the following problem

∂ₜ u + div f(u,p,t) = s(u,∇u,p,t)

f(u,t,i,j) where i is the line of f and j is the collumn of f

"""
function curve(s,i)
    ss = s
    if i==1
        return ss 
    elseif i==2
        return 10*cos(0.4*ss)
    else
        return 10*sin(0.4*ss)
    end
end

function curvep(s,i)
    return ForwardDiff.derivative(s->curve(s,i),s) 
end

function curvepp(s,i)
    return ForwardDiff.derivative(s->curvep(s,i),s) 
end

function normcurvepp(s)
    sqrt(sum(abs2,curvepp(s,i) for i in 1:3))
end


function curvepp_normal(s,i)
    normcpp = normcurvepp(s)
    return curvepp(s,i)/normcpp
end

function normcurvep(s)
    sqrt(sum(abs2,curvep(s,i) for i in 1:3))
end

function tang(s,i)
    normcp = normcurvep(s)
    return curvep(s,i)/normcp
end
function binormal(s,i)
    # tang ∧ curvepp_normal
    if i == 1
        a2 = tang(s,2)
        a3 = tang(s,3)
        b2 = curvepp_normal(s,2)
        b3 = curvepp_normal(s,3)
        return a2*b3-a3*b2
    elseif i==2
        a1 = tang(s,1)
        a3 = tang(s,3)
        b1 = curvepp_normal(s,1)
        b3 = curvepp_normal(s,3)
        return -a1*b3+a3*b1
    else
        a1 = tang(s,1)
        a2 = tang(s,2)
        b1 = curvepp_normal(s,1)
        b2 = curvepp_normal(s,2)
        return a1*b2-a2*b1
    end
end
function normal_curv(s,i)
    #binormal ∧ tang
    if i == 1
        a2 = binormal(s,2)
        a3 = binormal(s,3)
        b2 = tang(s,2)
        b3 = tang(s,3)
        return a2*b3-a3*b2
    elseif i==2
        a1 = binormal(s,1)
        a3 = binormal(s,3)
        b1 = tang(s,1)
        b3 = tang(s,3)
        return -a1*b3+a3*b1
    else
        a1 = binormal(s,1)
        a2 = binormal(s,2)
        b1 = tang(s,1)
        b2 = tang(s,2)
        return a1*b2-a2*b1
    end
end
function curvature(s)
    return normcurvepp(s)/normcurvep(s)
end
function er(θ,s,i)
    cos(θ)*normal_curv(s,i)+sin(θ)*binormal(s,i)
end
function eθ(θ,s,i)
    -sin(θ)*normal_curv(s,i)+cos(θ)*binormal(s,i)    
end

relu(x) = max(0,x)
@fastmath R0(p) = (p.y>11 || p.y<9) ? 2.0 : ( (p.x<-0.5 || p.x>0.5) ? 2.0 : 2.0+ 4*(p.y-9.0)*(p.y-11.0)*(p.x-0.5)*(p.x+0.5)/( (10-9.0)*(10-11.0)*(0-0.5)*(0+0.5)  ) )
function A0(p) 
    r0 = R0(p)
    return r0*r0/2.0
end
dA0θ(p) = ForwardDiff.derivative(θ->A0(Point(θ,p.y)),p.x)
dA0s(p) = ForwardDiff.derivative(s->A0(Point(p.x,s)),p.y)
P0(t) = t<0.0125 ? 13_332 : 0.0#sinpi(t/0.125)^2*2e4 : 0.0
@fastmath Ap(p,t) = (A0(p)*P0(t)/b(p) + sqrt(A0(p)))^2
@fastmath b(p) = 1e7*0.1/sqrt(2)
@fastmath P(A,p) = b(p)*(sqrt(A)-sqrt(A0(p)))/A0(p)
@fastmath Pp(A,p) = b(p)/(2*sqrt(A*A0(p)))
function f(u,p,t,i,j)
    a,qθ,qs = u
    A = a+A0(p)
    if i==1
        if j==1
            return qθ/A
        else
            return qs
        end
    elseif i==2
        if j==1
            return qθ^2/(2*A^2)+A*P(A,p)
        else
            return qθ*qs/A
        end
    else
        if j==1
            return qθ*qs/A^2
        else
            return qs^2/A-qθ^2/(2*A^2) + A*P(A,p)
        end
    end
end

function s(u,duθ,dus,p,t,i)
    a,qθ,qs = u
    A = a+A0(p)
    daθ = duθ[1]
    das = dus[1]
    dAθ = daθ+dA0θ(p)
    dAs = das+dA0s(p)
    C = curvature(p.y)
    R = sqrt(2*A)
    nu = 0.04
    if i==1
        return 0.0
    elseif i==2
        return 2/3*R*C*sin(p.x)*qs^2/A -2*11*nu*qθ/A + dAθ*P(A,p)
    else
        return -2/3*R*C*sin(p.x)*qθ*qs/A^2 -11*nu*qs/A + dAs*P(A,p)
    end
end

# we also need to define the eigen values of f called λ
@fastmath function λ(u,p,t,i,j)
    a,qθ,qs = u
    A = a+A0(p)
    if i==1
        if j==1
            return qθ/A^2
        else
            return qs/A
        end
    elseif i==2
        if j==1
            return sqrt(Pp(A,p))
        else
            return qs/A+sqrt(A*Pp(A,p))
        end
    else
        if j==1
            return -sqrt(Pp(A,p))
        else
            return qs/A-sqrt(A*Pp(A,p))
        end
    end
end

# About the initial condition now

function u0(p,i)
    if i==1
        return A0(p) - A0(p)
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
    return u[i]
end

function uright(u,p,t,i)
    return u[i]
end

function ubot(u,p,t,i)
    Ag,Qθg,Qsg = u
    Ag += A0(p) 
    A = Ap(p,t)
    if i==1
        return A-A0(p)
    elseif i==2
        return u[i]
    else
        return A*( Qsg/Ag + sqrt(Ag*Pp(Ag,p)) - sqrt(A*Pp(A,p)) )
    end
end

function utop(u,p,t,i)
    return u[i] # let's neumann everywhere here
end

#  here we will add a little function to get the cylinder back

function callback(df,cnt)
    # df is a dataframe (t,x,y,qθ/A,qs/A)
    dfcyl = @chain df begin
        @mutate(theta = x,s = y)
        @mutate(A = u1+!!A0(!!Point(theta,s)))
        @mutate(R = sqrt(2*A),P = P(A,!!Point(theta,s)))
        @mutate(ur = -dux2/(R*A)+2*u2/A^2*(dux1+dA0θ(!!Point(theta,s)))-duy3+u3/A*(duy1+dA0s(!!Point(theta,s))),
        utheta = 2*u2/A*1/R,
        us = u3/A
        )
        # now we should define the points
        @mutate(xc = !!curve(s,1) + R*!!er(theta,s,1),
        yc = !!curve(s,2) + R*!!er(theta,s,2),
        zc = !!curve(s,3) + R*!!er(theta,s,3))
        # we might also need the velocity vector
        @mutate(vx = 0*ur*!!er(theta,s,1) + utheta*!!eθ(theta,s,1) + us * !!tang(s,1),
        vy = 0*ur*!!er(theta,s,2) + utheta*!!eθ(theta,s,2) + us * !!tang(s,2),
        vz = 0*ur*!!er(theta,s,3) + utheta*!!eθ(theta,s,3) + us * !!tang(s,3))
        #now to simplify, the file should begin with (xc,yc,zc)
        @select(xc,yc,zc,vx,vy,vz,ur,utheta,us,Area=A,Radius=R,Pressure=P,theta,sax=s,time=t)
    end
    CSV.write("./datasCyl/data$cnt.csv",dfcyl)
    return nothing
end