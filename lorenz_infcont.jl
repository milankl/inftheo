using JLD
using PyPlot
using PyCall
@pyimport numpy as np

D = load("/home/kloewer/julia/lorenz_posit/dec_accuracy/data/lorenz_hr.jld")

x = D["xyz"][1,:]
y = D["xyz"][2,:]
z = D["xyz"][3,:]
##
N = length(x)

# histogram of x
bins = -20.1:0.3:20.1
C,bins = np.histogram(x,bins)
p = C/N         # pdf

# bin arrays
binsleft = bins[1:end-1]
binsright = bins[2:end]
binsmid = 1/2*(binsleft+binsright)

# conditional probability based on bit i
function conditional_histogram(x::Array,cond::Array,bins::Array,lag::Int)
    # cond is binary 0 or 1
    cond0 = (cond .== 0)[1:end-lag]
    cond1 = (cond .== 1)[1:end-lag]
    hist0,bins = np.histogram(x[lag+1:end][cond0],bins)
    hist1,bins = np.histogram(x[lag+1:end][cond1],bins)

    # normalize
    n0,n1 = sum(cond0),sum(cond1)
    p0,p1 = hist0/n0,hist1/n1
    q0,q1 = [n0,n1]/(n0+n1)     # a priori likelihood of bit being 0 or 1

    return p0,p1,q0,q1
end

function bits_signed_exp(x::Float32,i::Int)
    # converts the exponent bits from an unsigned integer to a signed integer
    if i == 1 || i >= 10    # the bit requested is not an exponent bit
        return bits(x)[i]
    else    # i is in [2,..,9]
        bias = 127
        k = parse(Int,"0b"*bits(x)[2:9])-bias    # signed exponent

        # subnormal numbers, i.e. x=0 means all exponent bits are 0
        if k == -bias   # i.e. x was 0
            k = 0
        end

        if i == 2   # sign of exponent
            if k < 0
                se = '1'
            else
                se = '0'
            end
        else
            return bits(abs(k))[end-9+i]    # the remaining 7 bits of the exponent
        end
    end
end

# entropy calculation
function entropy(p::Array)
    H = 0.
    for pi in p[:]
        if pi > 0.
            H -= pi*log2(pi)
        end
    end
    return H
end

function entropy_calc(p,p0,p1,q0,q1)

    Hx = entropy(p)
    Hb0 = entropy(p0)
    Hb1 = entropy(p1)

    # bitwise information content
    return Hx - q0*Hb0 - q1*Hb1
end

function information_content(x::Array{Float32,1},p::Array,bins::Array,lags::Array{Int})
    # preallocate
    nlags = length(lags)
    Icont = zeros(32,nlags)

    for (ilag,lag) in enumerate(lags)
        for bit in 1:32
            println("$ilag/$nlags,$bit/32")
            bi = [parse(Int,bits_signed_exp(xi,bit)) for xi in x]
            p0,p1,q0,q1 = conditional_histogram(x,bi,bins,lag)

            Icont[bit,ilag] = entropy_calc(p,p0,p1,q0,q1)
        end
    end
    return Icont
end


# preallocate
lags = Int.(round.(10.^(-1.5:0.1:2)./0.005))
Icont = information_content(x,p,bins,lags)

save("data/infcont_floats.jld","Icont",Icont)
