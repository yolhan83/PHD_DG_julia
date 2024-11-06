struct Point{T1,T2}
    x :: T1
    y :: T2
end
"""
    struct Element{T1<:Number,T2<:Number,T3<:Number,T4<:Number,T5<:Number,T6<:Number,T7<:Number,T8<:Number}
        center :: Point{T1,T2}
        h :: Tuple{T3,T4}
        id :: Int
        p ::Int
        gp ::SVector{Int,Point{T5,T6}}
        gw ::SVector{Int,T7}
        M ::SMatrix{T8}
    end
    call Element(center,h,id,p,gp,gw,M) or just Element(center,h,id,p)
"""
struct Element
    center :: Point{Float64,Float64}
    h :: Tuple{Float64,Float64}
    id :: Int
    p ::Int
    gp ::Vector{Point{Float64,Float64}}
    gw ::Vector{Float64}
    invM ::Matrix{Float64}
end
function Element(center,h,id,p) # add a constructor
    ndeg = (p+1)^2
    gp_dim,gw_dim = gauss(4*p+3)
    @tullio gp[i,j] := Point(gp_dim[i],gp_dim[j])
    @tullio gw[i,j] := gw_dim[i]*gw_dim[j]
    gp = gp[:]
    gw = gw[:]
    M = zeros(ndeg,ndeg)
    @tullio M[i,j] = gw[k]*phi(gp[k],i,p)*phi(gp[k],j,p)
    M .*= h[1]*h[2]/4
    invM = inv(M)
    return Element(center,h,id,p,gp,gw,invM)
end

"""
struct Face{T1<:Number,T2<:Number}
    left_id :: Int
    right_id :: Int
    normal ::SVector{2,T1}
    h ::T2
end

"""
struct Face
    el_left :: Element
    el_right :: Element
    normal ::Vector{Float64}
    center ::Point{Float64,Float64}
    h ::Float64
    gp ::Vector{Float64}
    gw ::Vector{Float64}
    tag ::Symbol
end
function Face(el_left,el_right,normal,center,h,deg,tag)
    gp,gw = gauss(4*deg+1+3)
    return Face(el_left,el_right,normal,center,h,gp,gw,tag)
end

"""
Make a simple mesh

"""

function make_mesh((Nx,Ny),deg,(ax,bx,ay,by);periodic=(false,false))
    h = ( (bx-ax)/Nx,(by-ay)/Ny )
    ele = Element[]
    fcs = Face[]
    for i in 1:Nx
        for j in 1:Ny
            center = Point( ax+(i-1)*h[1]+h[1]/2, ay+(j-1)*h[2]+h[2]/2 )
            id = (i-1)*Ny + j
            el = Element(center,h,id,deg)
            push!(ele,el)
            # take care of the bottom and left faces of elements
            if i>1
                normal = [1.0,0.0]                
                center = Point( ax+(i-1)*h[1], ay+(j-1)*h[2]+h[2]/2 )
                push!(fcs,Face(ele[id-Ny],el,normal,center,h[2],deg,:inner))
            else
                if !periodic[1]
                    normal = [1.0,0.0]                # fake left element
                    el_left = Element(Point(el.center.x-h[1],el.center.y),el.h,el.id,el.p )
                    center = Point( ax+(i-1)*h[1], ay+(j-1)*h[2]+h[2]/2 )
                    push!(fcs,Face(el_left,el,normal,center,h[2],deg,:left))
                end
            end
            if j>1
                    normal = [0.0,1.0]                
                    center = Point( ax+(i-1)*h[1]+h[1]/2, ay+(j-1)*h[2])
                    push!(fcs,Face(ele[id-1],el,normal,center,h[1],deg,:inner))
            else
                if !periodic[2]
                normal = [0.0,1.0]                # fake bottom element
                el_bot = Element(Point(el.center.x,el.center.y-h[2]),el.h,el.id,el.p )
                center = Point( ax+(i-1)*h[1]+h[1]/2, ay+(j-1)*h[2])
                push!(fcs,Face(el_bot,el,normal,center,h[1],deg,:bot))
                end
            end
            # take care of the top and right faces of elements (only need boundary ones)
            if i==Nx
                normal = [-1.0,0.0]
                if !periodic[1]
                    # fake right element
                    el_right = Element(Point(el.center.x+h[1],el.center.y),el.h,el.id,el.p )
                    center = Point( ax+(i-1)*h[1]+h[1], ay+(j-1)*h[2]+h[2]/2 )
                    push!(fcs,Face(el_right,el,normal,center,h[2],deg,:right))
                else
                    backid = Ny-(Nx*Ny-id)
                    el_right = ele[backid]
                    center = Point( ax+(i-1)*h[1]+h[1], ay+(j-1)*h[2]+h[2]/2 )
                    push!(fcs,Face(el_right,el,normal,center,h[2],deg,:inner))
                end
            end
            if j==Ny
                normal = [0.0,-1.0]
                if !periodic[2]
                    # fake top element
                    el_top = Element(Point(el.center.x,el.center.y+h[2]),el.h,el.id,el.p )
                    center = Point( ax+(i-1)*h[1]+h[1]/2, ay+(j-1)*h[2]+h[2] )
                    push!(fcs,Face(el_top,el,normal,center,h[1],deg,:top))
                else 
                    backid = id-Ny+1
                    el_top = ele[backid]
                    center = Point( ax+(i-1)*h[1]+h[1]/2, ay+(j-1)*h[2]+h[2] )
                    push!(fcs,Face(el_top,el,normal,center,h[1],deg,:inner))
                end
            end
        end
    end
    return ele,fcs
end

# define a mapping to [-1,1]x[-1,1]
function to_real(p,el::Element)
    # p is in [-1,1]x[-1,1]
    return Point(p.x*el.h[1]/2+el.center.x,p.y*el.h[2]/2+el.center.y)    
end
function to_map(p,el::Element)
    pt =  Point((p.x-el.center.x)*2/el.h[1], (p.y-el.center.y)*2/el.h[2])
    if pt.x >2 
        pt = Point(-1.0,pt.y)
    elseif pt.x < -2
        pt = Point(1.0,pt.y)
    end
    if pt.y >2
        pt = Point(pt.x,-1.0)
    elseif pt.y < -2
        pt = Point(pt.x,1.0)
    end
    return pt
end
function to_real(p,fc::Face)
    # here p is in [-1,1] and we need to recover the real point
    # we will suppose the unit teng vector is the direct rotation of the normal
    # now we remind we have the center of this face ! we need to moove from it in the teng direction by h/2*p
    return Point(fc.h/2*p*fc.normal[2] + fc.center.x,-fc.h/2*p*fc.normal[1] + fc.center.y)
end
function to_map(p,fc::Face)
    # here p is in the real and we need to recover a point in [-1,1]
    return ( (p.x-fc.center.x)*fc.normal[2] + (-fc.normal[1])*(p.y-fc.center.y) )*2/fc.h
end