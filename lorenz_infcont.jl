using JLD
using PyPlot
using StatsBase

D = load("/home/kloewer/julia/lorenz_posit/dec_accuracy/data/lorenz_hr.jld")

predictor = Float32.(D["xyz"][3,:])
predictand = Float32.(D["xyz"][2,:])

suffix = "zy"
println(suffix)

N = length(predictand)

# unconditional pdf, histogram of predictand
h = 0.3
#bins = -20.1:h:20.1     # for predictand x
bins = -27.3:h:27.3     # for predictand y
#bins = 0.0:h:50.0           # for predictand z
C = fit(Histogram,predictand,bins,closed=:left).weights
p = C/N         # pdf

println((p[1],p[end]))

# bin arrays
binsleft = bins[1:end-1]
binsright = bins[2:end]
binsmid = 1/2*(binsleft+binsright)

# conditional probability based on bit i
function conditional_histogram(x::Array,cond::Array,bins::StepRangeLen,lag::Int)
    # cond is binary 0 or 1
    cond0 = (cond .== 0)[1:end-lag]
    cond1 = (cond .== 1)[1:end-lag]
    hist0 = fit(Histogram,x[lag+1:end][cond0],bins,closed=:left).weights
    hist1 = fit(Histogram,x[lag+1:end][cond1],bins,closed=:left).weights

    # normalize
    n0,n1 = sum(cond0),sum(cond1)
    p0,p1 = hist0/n0,hist1/n1
    q0,q1 = [n0,n1]/(n0+n1)     # a priori likelihood of bit being 0 or 1

    return p0,p1,q0,q1
end

function bits_signed_exp(x::Float32)
    # converts the exponent bits from an unsigned integer to a signed integer
    all_bits = bits(x)  # unsigned exp

    bias = 127
    k = parse(Int,"0b"*bits(x)[2:9])-bias    # signed exponent

    # subnormal numbers, i.e. x=0 means all exponent bits are 0
    if k == -bias   # i.e. x was 0
        s_exp_bits = repeat("0",8)
    else
        if k < 0
            s0 = "1"
        else
            s0 = "0"
        end

        s_exp_bits = s0*bits(abs(k))[end-6:end]    # the remaining 7 bits of the exponent
    end
    return all_bits[1]*s_exp_bits*all_bits[10:end]
end

function entropy_calc(p,p0,p1,q0,q1)

    Hx = entropy(p,2)
    Hb0 = entropy(p0,2)
    Hb1 = entropy(p1,2)

    # bitwise information content
    return Hx - q0*Hb0 - q1*Hb1
end

function information_content(x::Array{Float32,1},y::Array{Float32,1},p::Array,bins::StepRangeLen,lags::Array{Int})
    # preallocate
    nlags = length(lags)
    Icont = zeros(32,nlags)

    # PREDICTOR X: convert only once to float32 bits
    B = zeros(Int8,32,length(x))
    for (xin,xi) in enumerate(x)
        for (ib,b) in enumerate(bits_signed_exp(xi))
            B[ib,xin] = parse(Int,b)
        end
    end

    for (ilag,lag) in enumerate(lags)
        for ibit = 1:32
            # use the predictand y here
            p0,p1,q0,q1 = conditional_histogram(y,B[ibit,:],bins,lag)
            Icont[ibit,ilag] = entropy_calc(p,p0,p1,q0,q1)
            println("$ilag/$nlags,$ibit/32")
        end
    end
    return Icont
end


# preallocate
lags = cat(1,[0,1,2,3,4,5],Int.(round.(10.^(-1.5:0.1:2)./0.005)))
Icont = information_content(predictor,predictand,p,bins,lags)

save("data/infcont_floats_$suffix.jld","Icont",Icont)
