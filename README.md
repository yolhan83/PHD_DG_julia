# PHD_DG_Julia

Welcome to my repository implementing the 2D Discontinous Galerkin Method (DGm) for advection system.

## How to run

To run the code simply follow those steps :

1. Install julia
2. git clone the repo using 
```bash 
git clone https://github.com/yolhan83/PHD_DG_julia/tree/main 
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

## How to make a test case

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

then you go to the "PHD_DG_julia.jl" file in "src" and change the import line 14 together with the parameter for the numerical simulation in main. Then, simply run as before.