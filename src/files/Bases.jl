@fastmath function phi(p,i,deg::Int) # special case when the Element is not created yet
    ligne = (i-1)÷(deg+1)
    col = (i-1)%(deg+1)
    res = one(p.x)
    @simd for j in 1:ligne
        res *= p.x
    end
    @simd for j in 1:col
        res *= p.y
    end
    return res
end

@fastmath function phi(p,i,el)
    ligne = (i-1)÷(el.p+1)
    col = (i-1)%(el.p+1)
    res = one(p.x)
    @simd for j in 1:ligne
        res *= p.x
    end
    @simd for j in 1:col
        res *= p.y
    end
    return res
end


@fastmath function phipx(p,i,el)
    ligne = (i-1)÷(el.p+1)
    if ligne == 0
        return zero(p.x)
    end
    col = (i-1)%(el.p+1)
    res = one(p.x)*ligne
    @simd for j in 1:ligne-1
        res *= p.x
    end
    @simd for j in 1:col
        res *= p.y
    end
    return res
end

@fastmath function phipy(p,i,el)
    col = (i-1)%(el.p+1)
    if col == 0
        return zero(p.x)
    end
    ligne = (i-1)÷(el.p+1)
    res = one(p.x)*col
    @simd for j in 1:ligne
        res *= p.x
    end
    @simd for j in 1:col-1
        res *= p.y
    end
    return res
end
