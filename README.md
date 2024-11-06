# PHD_DG_Julia

Welcome to my repository implementing the 2D Discontinous Galerkin Method (DGm) for advection system.

## How to run

To run the code simply follow those steps :

1. Install julia
2. git clone the repo using ```git clone https://github.com/yolhan83/PHD_DG_julia/tree/main ```
3. run ```julia -t auto --project setup.jl```
4. run ```julia -t auto --project -e "using PHD_DG_julia; PHD_DG_julia.main();"```

To clear the datas, 
run ```julia -t auto --project rmfiles.jl```