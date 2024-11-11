"""
make_df(U,param,t)
Make a dataframe from the solution at a certain time
"""
function make_df(U,param,t)
    ele = param[1]
    uncache = param[3] #@view(param[3][:,1])
    dunxcache = param[4] #@view(param[4][:,1])
    dunycache = param[5] #@view(param[5][:,1])
    ninc = size(U,3)
    numberofpoints = length(ele)*length(ele[1].gp)
    x = Vector{Float64}(undef,numberofpoints)
    y = Vector{Float64}(undef,numberofpoints)
    u = [ Vector{Float64}(undef,numberofpoints) for j in 1:ninc]
    dux = [ Vector{Float64}(undef,numberofpoints) for j in 1:ninc]
    duy = [ Vector{Float64}(undef,numberofpoints) for j in 1:ninc]
    Threads.@threads for n in eachindex(ele)
        el = ele[n]
        idx = Threads.threadid()
        uncacheid = @view(uncache[:,idx])
        dunxcacheid = @view(dunxcache[:,idx])
        dunycacheid = @view(dunycache[:,idx])
        for k in eachindex(el.gp,el.gw)
            p = el.gp[k]
            preal = to_real(p,el)
            cid = (n-1)*length(el.gp)+k
            x[cid] = preal.x
            y[cid] = preal.y
            eval_u!(uncacheid,U,p,el)
            eval_dux!(dunxcacheid,U,p,el)
            eval_duy!(dunycacheid,U,p,el)
            for j in 1:ninc
                u[j][cid] = uncacheid[j]
                dux[j][cid] = dunxcacheid[j]
                duy[j][cid] = dunycacheid[j]
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
        @inbounds  @simd for j in axes(U,3)
            eval_u!(@view(un[:,tid]),U,Point(0.0,0.0),el)
            @inbounds  @simd for d in 1:2
                m[tid] = max(m[tid],abs( λ(@view(un[:,tid]),p,t,j,d) ))
            end
        end
    end
    global_speed = maximum(m)
    return cfl*hmin/(2*degmax+1)/global_speed , global_speed
end