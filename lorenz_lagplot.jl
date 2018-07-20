using JLD
using PyPlot
using StatsBase

#D = load("/home/kloewer/julia/lorenz_posit/dec_accuracy/data/lorenz_hr.jld")

path = "/local/kloewer/julsdata/"
runs = [5]

# LOAD DATA
ncu = NetCDF.open(path*"run"*@sprintf("%04d",runs[1])*"/u.nc")
predictor = reshape(ncu.vars["u"][5,50,:],40019)
predictand = reshape(ncu.vars["u"][5,50,:],40019)

#x = Float32.(D["xyz"][1,:])
#y = Float32.(D["xyz"][2,:])
#z = Float32.(D["xyz"][3,:])
#
N = length(predictand)

# histogram of x
#bins = -20.1:0.3:20.1
h = 0.05
#bins = 0:0.3:50
bins = -1.5:h:1.5
C = fit(Histogram,predictand,bins,closed=:left).weights
p = C/N         # pdf

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

#

#using SigmoidNumbers
#P = Posit{32,2}

ion()
fig,axs = subplots(1,3,sharex=true,sharey=false,figsize=(10,3))

i = 1
for (ax,lag) in zip([axs[1],axs[2],axs[3]],[0,8,40])

    bi = [parse(Int,bits_signed_exp(xi,i)) for xi in predictor]
    #bi = [parse(Int,bits(P(xi))[i]) for xi in predictor]

    p0,p1,q0,q1 = conditional_histogram(predictand,bi,bins,lag)

    Hx = entropy(p,2)
    Hb0 = entropy(p0,2)
    Hb1 = entropy(p1,2)

    # bitwise information content
    Ib = Hx - q0*Hb0 - q1*Hb1

    # plotting
    function f2(x)
        return @sprintf("%.3f",x)
    end

    #
    ax[:plot](binsmid,p,drawstyle="steps-post",label=L"$p(z)$")
    ax[:plot](binsmid,p0,drawstyle="steps-post",label=L"$p(z|x_3(t-τ) = 0)$")
    ax[:plot](binsmid,p1,drawstyle="steps-post",label=L"$p(z|x_3(t-τ) = 1)$")
    ax[:text](0.7,0.85,"I = $(f2(Ib))",transform=ax[:transAxes])
end

axs[1][:set_title]("instantaneous (τ = 0)")
axs[2][:set_title]("short forecast (τ = 2d)")
axs[3][:set_title]("long forecast (τ = 10d)")

axs[1][:set_xlabel](L"$z$")
axs[2][:set_xlabel](L"$z$")
axs[3][:set_xlabel](L"$z$")
axs[1][:legend](loc=4,fontsize=7)

axs[1][:set_title]("a",loc="left",fontweight="bold")
axs[2][:set_title]("b",loc="left",fontweight="bold")
axs[3][:set_title]("c",loc="left",fontweight="bold")


tight_layout()
#savefig("figs/inf_cont_xz_P_$i.pdf")
#savefig("bitwise_inf.pdf")

#close(fig)
