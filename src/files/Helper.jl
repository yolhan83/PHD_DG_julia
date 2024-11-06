"""
make_df(U,param,t)
Make a dataframe from the solution at a certain time
"""
function make_df(U,param,t)
    ele = param[1]
    uncache = @view(param[3][:,1])
    dunxcache = @view(param[4][:,1])
    dunycache = @view(param[5][:,1])
    x = Float64[]
    y = Float64[]
    u = [ Float64[] for j in 1:ninc]
    dux = [ Float64[] for j in 1:ninc]
    duy = [ Float64[] for j in 1:ninc]
    for el in ele
        id = el.id
        for k in eachindex(el.gp,el.gw)
            p = el.gp[k]
            preal = to_real(p,el)
            push!(x,preal.x)
            push!(y,preal.y)
            eval_u!(uncache,U,p,el)
            eval_dux!(dunxcache,U,p,el)
            eval_duy!(dunycache,U,p,el)
            for j in 1:ninc
                push!(u[j],uncache[j])
                push!(dux[j],dunxcache[j])
                push!(duy[j],dunycache[j])
            end
        end
    end
    return DataFrame(:x=>x,:y=>y,[Pair(Symbol("u$j"),u[j]) for j in 1:ninc]..., [Pair(Symbol("dux$j"),dux[j]) for j in 1:ninc]...,[Pair(Symbol("duy$j"),duy[j]) for j in 1:ninc]... ,:t => t.*ones(length(x)))
end

"""
adapt_dt(U,param,t;cfl = 0.9)
Adapt the time step using the CFL condition for RKDG method
"""
function adapt_dt(U,param,t;cfl = 0.9) 
    ele = param[1]
    hmin = minimum(el.h[1]*el.h[2]/(2*el.h[1]+2*el.h[2]) for el in ele)
    degmax = maximum(el.p for el in ele)
    un = param[3]
    m = zeros(size(un,2))
    Threads.@threads for el in ele
        id = el.id
        tid = Threads.threadid()
        p = el.center
        for j in axes(U,3)
            eval_u!(@view(un[:,tid]),U,Point(0.0,0.0),el)
            for d in 1:2
                m[tid] = max(m[tid],abs( λ(@view(un[:,tid]),p,t,j,d) ))
            end
        end
    end
    global_speed = maximum(m)
    return cfl*hmin/(2*degmax+1)/global_speed , global_speed
end