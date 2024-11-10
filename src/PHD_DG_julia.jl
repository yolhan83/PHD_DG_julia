module PHD_DG_julia

try
    using Tullio,QuadGK,DataFrames,ForwardDiff,CSV,TidierData,Plots,ThreadsX 
catch
    throw(ErrorException("run 'julia -t auto --project setup.jl' first"))
end
include("files/Bases.jl") # for basis functions
include("files/mesh.jl") # for mesh (really bad may use a librairy one day)
include("files/Helper.jl")
include("files/RKsolver.jl") # for RK solver
include("files/DG.jl") # for DG semi-solver

include("TCs/BF.jl") # for advection model

function main()
    include("./src/NumParam.jl")
    param = get_param(N...,deg,ax,bx,ay,by;periodic = (true,false)) 
    prob = get_prob(param,ninc,deg,tspan) 
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
    
        @time solve(prob, method;nsave=30,callback = callback,verbose = true); 
    end    
end
end # module PHD_DG_julia
