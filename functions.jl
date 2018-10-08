function continuous_idx(v::Array,x::Real)
    if v != sort(v)
        throw(error("Only monotonically increasing arrays are allowed."))
    end

    i1 = find_closest(v,x)
    if v[i1] > x
        i1 -= 1
    end
    return i1+(x-v[i1])/(v[i1+1]-v[i1])
end

function find_closest(v::Array,x::Real)
    return argmin(abs.(v .- x))
end
