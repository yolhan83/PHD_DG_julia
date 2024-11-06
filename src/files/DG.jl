"""
get_U0(ele,ninc)
project a function u0(x) on the elements
"""
function get_U0(ele,ninc)
    N = length(ele)
    deg = maximum(el.p for el in ele)
    U0 = zeros((deg+1)^2,N,ninc)
    V = zeros((deg+1)^2)
    for el in ele
        id = el.id
        for j in axes(U0,3)
            V .= 0.0
            for k in eachindex(el.gp,el.gw)
                p = el.gp[k]
                preal = to_real(p,el)
                for i in axes(U0,1)
                    V[i] += el.h[1]*el.h[2]/4*el.gw[k]*u0(preal,j)*phi(p,i,el)
                end
            end
            @tullio U0[i1,$id,$j] = el.invM[i1,i2]*V[i2]
        end
    end
    return U0
end

"""
get_param(Nx,Ny,deg,ax,bx,ay,by)
make the mesh and some caching for our DG method, store in parameters of the later ODEprob
"""
function get_param(Nx,Ny,deg,ax,bx,ay,by)
    ele,fcs = make_mesh((Nx,Ny),deg,(ax,bx,ay,by));
    nth = Threads.nthreads()
    return (ele,fcs,zeros(3,nth),zeros(3,nth),zeros(3,nth),zeros((deg+1)^2));
end

"""
get_prob(param,ninc,deg,tspan)
make the ODEprob
"""
function get_prob(param,ninc,deg,tspan)
    ele = param[1]
    U0 = get_U0(ele,ninc)
    prob = ODEprob(ode!,U0,tspan,param)
end

"""
eval_u!(u,U,p,el)
get the value of u at map point p ∈ [-1,1]^2 in element el and store in u
"""
function eval_u!(u,U,p,el)
    # get the value of u at map point p ∈ [-1,1]^2 in element el
    id = el.id
    u .= 0.0
    @inbounds @simd for j in eachindex(u)
        @inbounds @simd for i in axes(U,1)
            u[j] += U[i,id,j]*phi(p,i,el)
        end
    end 
end

"""
eval_dux!(dux,U,p,el)
get the value of du/dx at map point p ∈ [-1,1]^2 in element el and store in dux
"""
function eval_dux!(dux,U,p,el)
    id = el.id
    dux .= 0.0
    @inbounds @simd for j in eachindex(dux)
        @inbounds @simd for i in axes(U,1)
            dux[j] += 2/el.h[1]*U[i,id,j]*phipx(p,i,el)
        end
    end
end

"""
eval_duy!(duy,U,p,el)
get the value of du/dy at map point p ∈ [-1,1]^2 in element el and store in duy
"""
function eval_duy!(duy,U,p,el)
    id = el.id
    duy .= 0.0
    @inbounds @simd for j in eachindex(duy)
        @inbounds @simd for i in axes(U,1)
            duy[j] += 2/el.h[2]*U[i,id,j]*phipy(p,i,el)
        end
    end
end

"""
volume_integral!(dU,U,param,t)
Calculate the volume integral part of the DG method and store in dU
"""
function volume_integral!(dU,U,param,t)
    ele = param[1]
    un_cache = param[3]
    dunx = param[4]
    duny = param[5]
    Threads.@threads for el in ele
        # take the id
        id = el.id
        tid = Threads.threadid()
        # let's work on a specific gauss point
        @inbounds for k in eachindex(el.gp,el.gw) # again make sure we're inbound and optimize the integration
            # we will need the current value of u at the gauss point
            p = el.gp[k]
            preal = to_real(p,el)
            eval_u!(@view(un_cache[:,tid]),U,p,el)
            eval_dux!(@view(dunx[:,tid]),U,p,el)
            eval_duy!(@view(duny[:,tid]),U,p,el)
            # we compute the volume integral ∫ f ⋅ ∇ϕ
            @inbounds for j in axes(dU,3)
                F1 = f(@view(un_cache[:,tid]),preal,t,j,1)
                F2 = f(@view(un_cache[:,tid]),preal,t,j,2)
                S = s(@view(un_cache[:,tid]),@view(dunx[:,tid]),@view(duny[:,tid]),preal,t,j)
                @inbounds for i in axes(dU,1)
                    dU[i,id,j] += el.gw[k]*(el.h[2]/2*F1*phipx(p,i,el) + el.h[1]/2*F2*phipy(p,i,el))
                    dU[i,id,j] += el.gw[k]*(el.h[1]*el.h[2]/4*S*phi(p,i,el))
                end
            end
        end
    end
end

"""
Rusanov_flux(uleft,uright,normal,p,t,j)
Calculate the Rusanov flux of the jth variable at a map point on the interface
"""
@fastmath function Rusanov_flux(uleft,uright,normal,p,t,j,c)
    mean_flux = (f(uleft,p,t,j,1)+f(uright,p,t,j,1))*normal[1]/2 + (f(uleft,p,t,j,2) + f(uright,p,t,j,2))*normal[2]/2
    jump = uright[j] - uleft[j]
    local_speed = c # maximum(abs( λ(uleft,p,t,jj,1)*normal[1] + λ(uleft,p,t,jj,2)*normal[2] ) for jj in eachindex(uleft))
    # @info c
    return mean_flux - local_speed*jump
end

"""
surface_integral!(dU,U,param,t)
Calculate the surface integral part of the DG method and store in dU
"""
function surface_integral!(dU,U,param,t,c)
    ele = param[1]
    fcs = param[2]
    unleft = param[3]#@view(param[3][:,1])
    unright = param[4]#@view(param[4][:,1])
    uncache = param[5]#@view(param[5][:,1])
    for fc in fcs
        tid = Threads.threadid()
        unleftid= @view(unleft[:,tid])
        unrightid= @view(unright[:,tid])
        uncacheid= @view(uncache[:,tid])
        el_left = fc.el_left
        idleft = el_left.id
        el_right = fc.el_right
        idright = el_right.id
        normal = fc.normal
        # let's work on a specific gauss point
        for k in eachindex(fc.gp)
            p = fc.gp[k]
            # we first need the real point
            preal = to_real(p,fc)
            # we also need the ref point in the right/left element
            pright = to_map(preal,el_right)
            pleft = to_map(preal,el_left)
            # now the values of u at the gauss point, here we need the tag, also note that the right element always exists
            eval_u!(unrightid,U,pright,el_right)
            if fc.tag == :left
                # this means we are enterring the domain at left
                eval_u!(uncacheid,U,pleft,el_right)
                for j in eachindex(unleftid)
                    unleftid[j] = uleft(uncacheid,preal,t,j) 
                end
            elseif fc.tag == :right
                # this means we are leaving the domain at right
                eval_u!(uncacheid,U,pleft,el_right)
                for j in eachindex(unleftid)
                    unleftid[j] = uright(uncacheid,preal,t,j) 
                end
            elseif fc.tag == :top
                # this means we are leaving the domain at top
                eval_u!(uncacheid,U,pleft,el_right)
                for j in eachindex(unleftid)
                    unleftid[j] = utop(uncacheid,preal,t,j) 
                end
            elseif fc.tag == :bot
                # this means we are entering the domain at bottom
                eval_u!(uncacheid,U,pleft,el_right)
                for j in eachindex(unleftid)
                    unleftid[j] = ubot(uncacheid,preal,t,j) 
                end
            elseif fc.tag == :inner
                eval_u!(unleftid,U,pleft,el_left)
            else
                error("unknown tag")
            end
            # ok we have everything to calculate the numerical flux !
            for j in axes(dU,3)
                flux = Rusanov_flux(unleft,unright,normal,preal,t,j,c)
                for i in axes(dU,1)
                    if fc.tag == :inner
                        @inbounds dU[i,idleft,j] -= fc.h/2*fc.gw[k]*flux*phi(pleft,i,fc.el_left)
                    end
                    @inbounds dU[i,idright,j] += fc.h/2*fc.gw[k]*flux*phi(pright,i,fc.el_right)
                end
            end
        end
    end
    nothing
end


"""
ode!(dU,U,param,t)
Calculate the semi-dicretization of the ODE as U'(t) = ode(U(t),p,t) and store in dU.
"""
function ode!(dU,U,param,t,c)
    dU .= 0.0 # init at zeros
    volume_integral!(dU,U,param,t) # we compute the volume integral
    surface_integral!(dU,U,param,t,c) # we compute the surface integral
    ele = param[1]
    DOF_cache = param[6]
    for el in ele
        id = el.id
        invM = el.invM
        for j in axes(dU,3)
            @tullio DOF_cache[i1] = invM[i1,i2]*dU[i2,$id,$j]
            @views dU[:,id,j] .= DOF_cache
        end
    end
    return nothing
end


