module PHD_DG_julia


using Tullio,QuadGK,DataFrames,ForwardDiff,CSV,TidierData,Plots,ThreadsX
include("src/files/Bases.jl") # for basis functions
include("src/files/mesh.jl") # for mesh (really bad may use a librairy one day)
include("src/files/Helper.jl")
include("src/files/RKsolver.jl") # for RK solver
include("src/files/DG.jl") # for DG semi-solver

include("./TCs/BF.jl") # for advection model


N = (32,64)
deg=2
ax,bx = -pi,pi
ay,by = 0.0,20.0
param = get_param(N...,deg,ax,bx,ay,by) 


ninc = 3
tspan = (0.0,0.04)
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

du = similar(prob.u0)
@profview for i in 1:10 ode!(du, prob.u0, param, 0.0,c) end

dir = "./datasCyl"
files = readdir(dir)
using ThreadsX
df = reduce(vcat,ThreadsX.map(file->CSV.read(joinpath(dir,file),DataFrame),files))
names(df)
df_1d_time = @chain df begin
    @group_by(time,sax)
    @summarise(across(1:15,mean))
    @arrange(time,sax)
end

using Plots

@gif for (i,df) in enumerate(df_1d_time)
    plot(df.sax,df.us_mean,legend=false)
end

end # module PHD_DG_julia
