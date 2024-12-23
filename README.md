# PHD_DG_Julia

Welcome to my repository implementing the 2D Discontinous Galerkin Method (DGm) for advection system.

## How to run ?

To run the code simply follow those steps :

1. Install julia
2. git clone the repo using 
```bash 
git clone https://github.com/yolhan83/PHD_DG_julia/
```
3. cd the direcrory
```bash 
cd PHD_DG_julia
```
3. setup the project (only ones) 
```bash
julia -t auto --project setup.jl
```
4. Create the folder "datas" and the folder "datasCyl" (only for blood flow)
5. run 
```bash 
julia -t auto --project -e "using PHD_DG_julia; PHD_DG_julia.main();"
```

To clear the datas, 
run 
```bash 
julia -t auto --project rmfiles.jl
```

## How to make a test case ?

A test case need to implement the follwing functions

```julia
f(u::AbstractVector,p::Point,t::Number,i::Int,j::Int) ::Number

s(u::AbstractVector,dux::AbstractVector,duy::AbstractVector,p::Point,t::Number,i::Int) ::Number

λ(u::AbstractVector,p::Point,t::Number,i::Int,j::Int) ::Number

u0(p::Point,i::Int) ::Number

uleft(u::AbstractVector,p::Point,t::Number,i::Int) ::Number

uright(u::AbstractVector,p::Point,t::Number,i::Int) ::Number

ubot(u::AbstractVector,p::Point,t::Number,i::Int) ::Number

utop(u::AbstractVector,p::Point,t::Number,i::Int) ::Number

callback(df::DataFrame,cnt::Int) ::Nothing
```

this package will then solve the following problem,
```math
\begin{cases}
    \partial_t u + div f(u,p,t) &=& s(u,\nabla u,p,t)\\
    u(p,t=0) &=& u_0(p)\\
    u((ax,y),t) &=& u_{left}(t)\\
    u((bx,y),t) &=& u_{right}(t)\\
    u((x,ay),t) &=& u_{bot}(t)\\
    u((x,by),t) &=& u_{top}(t)
\end{cases} \quad  ax \le x \le bx,\, ay \le y \le by, \, 0 \le t \le T
```
with $u\in \mathbb{R}^n$, $f\in \mathbb{R}^{n\times d}$, $s\in \mathbb{R}^n$. Also, $\lambda \in \mathbb{R}^{n\times d}$ is the eigen values of $\nabla_u f$. Note that $n$ is the number of unknown and the number of equations too here and $d$ is the dimension $d=2$ evereywhere here.


Finally, you go to the "PHD_DG_julia.jl" file in "src" and change the import line 14 together with the parameter for the numerical simulation in the NumParam file in src. Then, simply run as before.

## How to view my datas ?

You can view your datas using Visit and Paraview or directly in julia.

For visit : open the datas folder where it was saved and open them using PlainText, check the box "first row as variable names" and set x and y coordinates as 0 and 1 index respectivly. 

For paraview : open the datas forlder as csv view, use filter "datas to point" and configure x,y to be the coordinates.

For julia : in a julia repl, you will have to install DataFrames, CSV and Plots using
```julia
using Pkg
Pkg.add(["Plots","DataFrames","CSV"])
```
Then run the follwing code,
```julia
using Plots,DataFrames,CSV
dir = "./datas" # or the complete path
files = readdir(dir)
df = reduce(vcat,map(f -> CSV.read(f,DataFrame),files))
heatmap(df.x,df.u,df.u1)
```

Adapt this to datasCyl where (xc,yc,zc) are the coordinates (index 0,1,2).

## Remarks

1. The boundary conditions takes $u$ as an argument, this $u$ is the one coming from the trace of the solution at the face. This allows for simple neumann conditions.
2. Boundary conditions are all weakly impose, meaning, if they're not "physical" they won't be forced. User should then look for "hyperbolicly right" boundary conditions. If you think they should be right for a specific case, I will be glad to talk about it !

## Future plan

This should begin after the phd is fully over, also, tests should be added, and, a seperate package for the features should be made so that this one would just be an ease to run for julia non initiated.

1. Support any Meshes.jl mesh or T8Code.jl mesh.
2. After (1) is solved, implement what's needed for MPI.jl to be used, a first step would be to allow for some kind of distributive programing
3. Change to use DifferentialEquations.
4. After (3), implement a parabolic term with Additive methods.

## Notes

All those functionlity with much more (any P4Est, T8Code mesh, AMR callbacks at any point (RK substep, step, save), adaptive mesh, cpu, mpi, gpu .... ) are implemented within Trixi.jl, use that is you want DG at full power.