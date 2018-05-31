using PyPlot
using JLD
using PyCall
@pyimport cmocean.cm as cmaps

cd("/home/kloewer/julia/inftheo/")

Icont = load("data/infcont_floats.jld")["Icont"]

Icont[Icont.<=0] = NaN
lIcont = log.(Icont)
lags = Int.(round.(10.^(-1.5:0.1:2)./0.005))

bitvec = repmat(-.5 + (1:size(Icont)[1]),1,size(Icont)[2])
lagvec = repmat(1:size(Icont)[2],1,size(Icont)[1])'

ylabelsign = [L"sign"]
ylabelexp = [LaTeXString('$'*"e_$i"*'$') for i in 1:8]
ylabelfrac = [LaTeXString('$'*"f_{$i}"*'$') for i in 1:23]

ylabels = cat(1,ylabelsign,ylabelexp,ylabelfrac)

fig,(ax1,ax2) = subplots(1,2,figsize=(10,5),sharey=true)
Q1 = ax1[:pcolormesh](lagvec,bitvec,lIcont,vmin=-10,vmax=1,cmap="cubehelix_r")
Q2 = ax2[:pcolormesh](lagvec,bitvec,lIcont,vmin=-10,vmax=1,cmap="cubehelix_r")

ax1[:invert_yaxis]()
ax1[:set_yticks](Array(1:size(Icont)[1]))
ax1[:set_yticklabels](ylabels)
colorbar(ax=(ax1,ax2),Q2)


xtiks = [find(lags*0.005 .== x) for x in [0.1,1.,10.,100.]]
ax1[:set_xticks](xtiks)
ax1[:set_xticklabels]([L"10^{-1}",L"10^{0}",L"10^{1}",L"10^{2}"])
ax1[:set_xlabel]("forecast time (mtu)")

ax1[:plot]([0,36],[1.5,1.5],"k",lw=2)
ax1[:plot]([0,36],[9.5,9.5],"k",lw=2)

ax2[:plot]([0,36],[1.5,1.5],"k",lw=2)
ax2[:plot]([0,36],[9.5,9.5],"k",lw=2)

ax1[:set_ylim](32.5,.5)
ax1[:set_xlim](1,36)
ax2[:set_xlim](1,36)
ax1[:set_title]("Bit-wise information content, 32-bit Floats")
