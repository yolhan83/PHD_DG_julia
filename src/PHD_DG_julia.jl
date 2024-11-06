module PHD_DG_julia

try
    using Tullio,QuadGK,DataFrames,ForwardDiff,CSV,TidierData,Plots,ThreadsX 
catch
    throw(ErrorException("run 'julia -t auto --project setup.jl first'"))
end
include("files/Bases.jl") # for basis functions
include("files/mesh.jl") # for mesh (really bad may use a librairy one day)
include("files/Helper.jl")
include("files/RKsolver.jl") # for RK solver
include("files/DG.jl") # for DG semi-solver

include("TCs/BF.jl") # for advection model

function main()
    N = (16,32)
    deg=2
    ax,bx = -pi,pi
    ay,by = 0.0,20.0
    param = get_param(N...,deg,ax,bx,ay,by) 
    
    ninc = 3
    tspan = (0.0,0.01)
    prob = get_prob(param,ninc,deg,tspan) 
    
    time_order = deg+1
    method = RKMethod(time_order)()
    dtmax,c = adapt_dt(prob.u0,param,0.0) 
    
    begin
        dir = "./datas"
        file = readdir(dir)
        for f in file
            rm(joinpath(dir,f))
        end
    
        dir = "./datasCyl"
        file = readdir(dir)
        for f in file
            rm(joinpath(dir,f))
        end
    
        @time solve(prob, method;nsave=30,callback = cylinder,verbose = true); 
    end    
end
end # module PHD_DG_julia
