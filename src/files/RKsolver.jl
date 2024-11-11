"""
struct RKmethod{T<:Real}  # all struct should be typed in julia 
    s :: Int # number of stages
    a :: Matrix{T} # matrix of coefficients
    b :: Vector{T} # vector of coefficients
    c :: Vector{T} # vector of coefficients
    b̂ :: Vector{T} # vector of coefficients 
end
Stores a RK method
"""
struct RKmethod{T<:Real}  # all struct should be typed in julia 
    s :: Int # number of stages
    a :: Matrix{T} # matrix of coefficients
    b :: Vector{T} # vector of coefficients
    c :: Vector{T} # vector of coefficients
end
# create the one you want
"""
RK1()
_____
0 | 0
_____
    1
"""
function RK1()
    a = [0;;]
    b = [1]	
    c = [0]
    RKmethod(1, a, b, c)
end
"""
RK2()
_____
0 | 0
1 | 1 0
_____
  | 1/2 1/2
"""
function RK2()
    a = [0 0
        1//1 0]
    b = [1//2, 1//2]
    c = [0//1,1]
    RKmethod(2, a, b, c)
end
"""
RK3()
________________
0   | 0
1   | 1   0
1/2 | 1/4 1/4 0
________________
    | 1/6 1/6 2/3
"""
function RK3()
    a = [0 0 0
        1 0 0
        1//4 1//4 0]
    b = [1//6, 1//6,2//3]
    c = [0,1,1//2]
    RKmethod(3, a, b, c)
end
"""
RK4()
_____________________
0   | 0
1/2 | 1/2 0
1/2 | 0   1/2 0
1   | 0   0   1   0
_____________________
    | 1/6 1/3 1/3 1/6
"""
function RK4()
    a = [0 0 0 0
        1//2 0 0 0
        0 1//2 0 0
        0 0 1 0]
    b = [1//6, 1//3, 1//3, 1//6]
    c = [0, 1//2, 1//2, 1]
    RKmethod(4, a, b, c)
end
"""
RK5()
"""
function RK5()
    a = [0 0 0 0 0 0 0
    1//5 0 0 0 0 0 0
    3//40 9//40 0 0 0 0 0
    44//45 -56//15 32//9 0 0 0 0
    19372//6561 -25360//2187 64448//6561 -212//729 0 0 0
    9017//3168 -355//33 46732//5247 49//176 -5103//18656 0 0
    35//384 0 500//1113 125//192 -2187//6784 11//84 0]
    b = [35//384,0,500//1113,125//192,-2187//6784,11//84,0]
    c = [0,1//5,3//10,4//5,8//9,1,1]
    RKmethod(7, a, b, c)
end
"""
RKMethod(time_order::Int)
Choose the right method
"""
function RKMethod(time_order::Int)
	if time_order==1
		method = RK1
	elseif time_order==2
		method = RK2
	elseif time_order==3
		method = RK3
	elseif time_order==4
		method = RK4
	else
		method = RK5
	end
    return method
end

# store everything about your problem
"""
struct ODEprob{F<:Function,T1<:AbstractArray,T2,T3}
    f::F
    u0 ::T1
    tspan ::T2
    param ::T3 
end
Store an ODE problem u'(t) = f(u(t),param,t), u(t0) = u0
"""
struct ODEprob{F<:Function,T1<:AbstractArray,T2,T3}
    f::F
    u0 ::T1
    tspan ::T2
    param ::T3 
end
"""
step!(ode!, u, cache, t, dt,param, method)
Make a step with the RK method choosen, note that the cache should be s+1 vector of same size as u.
"""
function step!(ode!, y, cache, t, dt,param, method,cc)
    # load the method parameters
    s = method.s
    a = method.a
    b = method.b
    c = method.c
    K = @view(cache[1:end-1])
    ysub = cache[end]
    # main loop cache[1:s] are the substeps and cache[end] is an helpful cache
    for i in eachindex(K,c)
        # calculate y + dt*∑( a_{ij}K_j ) and store it in cache[end]
        for k in eachindex(ysub,y,K[i])
            ysub[k] = y[k]
            for j in 1:i-1
                ysub[k] += dt*a[i,j] * K[j][k]
            end
        end
        # cache[i] = K_i = f(ysub,t+dt*c_i)
        ode!(K[i], ysub, param, t+c[i]*dt,cc)
    end 
    # deduce your next step
    for i in 1:s
        for k in eachindex(y)
            y[k] += dt * b[i] * cache[i][k]
        end
    end
    nothing
end
"""
solve(ODEprob, method, dtmax;nsave=100,callback = nothing;verbose = false)
Solve an ODE problem using the method choosen. The solution will be stored nsave times and the callback used After every save point. verbose will give time info at each step
"""
function solve(ODEprob, method;nsave=100,callback = nothing,verbose = false)
    # load problem
    ode! = ODEprob.f
    tspan = ODEprob.tspan
    y0 = ODEprob.u0
    param = ODEprob.param
    # about time
    t = tspan[1]
    dt,c = adapt_dt(y0,param,0.0) 
    saveat = (tspan[2]-tspan[1])/(nsave-1)
    # about the solution
    y = copy(y0) # avoid cahing the yo in ODEprob
    # cache reduce allocation
    cache = [ similar(y0) for _ in 1:method.s+1]
    ode!(cache[1], y0, param, t,c)
    time_fun = time_ns()
    ode!(cache[1], y0, param, t,c)
    time_fun = time_ns() - time_fun
    estimated_time = time_fun/1e9*(tspan[2]-tspan[1])/dt*method.s/60
    @info "estimated time : $(floor(Int,estimated_time)) minutes and $(floor(Int,(estimated_time-floor(Int,estimated_time))*60)) seconds"
    # main loop
    cnt = 1
    df = make_df(y,param,t)
    CSV.write("./datas/data$cnt.csv",df)
    if !isnothing(callback)
        callback(df,cnt)
    end
    while t <= tspan[2]
        # try
            step!(ode!,y,cache,t,dt,param,method,c)
        #     catch e
        #         @error e
        #         break
        # end
        t += dt
        if verbose
            println("time : $t")
        end
        if t > cnt*saveat
            cnt += 1
            df = make_df(y,param,t)
            CSV.write("./datas/data$cnt.csv",df)
            if !isnothing(callback)
                callback(df,cnt)
            end
            dt,c = adapt_dt(y,param,t) 
            @info cnt
        end
    end
    return nothing
end