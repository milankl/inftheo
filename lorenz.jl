using JLD
using PyCall
@pyimport numpy as np

D = load("/home/kloewer/julia/lorenz_posit/dec_accuracy/data/lorenz_hr.jld")

x = D["xyz"][1,:]
y = D["xyz"][2,:]
z = D["xyz"][3,:]

N = length(x)

function entropy(p::Array)
    H = 0.
    for pi in p[:]
        if pi > 0.
            H -= pi*log2(pi)
        end
    end
    return H
end

# histogram of x and entropy of x
bins = -20.1:0.3:20.1
C,bins = np.histogram(x,bins)

binsleft = bins[1:end-1]
binsright = bins[2:end]
binsmid = 1/2*(binsleft+binsright)

pplus = sum(sign.(x) .== 1)/N
pminus = 1-pplus

p = C/N
cplus = p/pplus
cminus = p/pminus

# left edge defines the sign within bin - careful with bin edges!
cplus[sign.(binsmid) .== -1] = 0
cminus[sign.(binsmid) .== 1] = 0

Hx = entropy(p)
Hxplus = entropy(cplus)
Hxminus = entropy(cminus)

plot(binsmid,p)
plot(binsmid,cplus)
plot(binsmid,cminus)
title(Hx)
