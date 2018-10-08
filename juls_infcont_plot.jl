using PyPlot
using JLD2
using PyCall
using FileIO
@pyimport cmocean.cm as cmaps

cd("/home/kloewer/julia/inftheo/")

file1 = "data/juls/infcont_floats_uu.jld2"
file2 = "data/juls/infcont_floats_uu4.jld2"

Icont1 = load(file1)["Icont"]
Icont2 = load(file2)["Icont"]

# remove nans and convert to base log10
Icont1[Icont1.<=0] .= NaN
lIcont1 = log10.(Icont1)

Icont2[Icont2.<=0] .= NaN
lIcont2 = log10.(Icont2)

lags1 = load(file1)["lags"]
lags2 = load(file2)["lags"]

dt1 = load(file1)["dt"]
dt2 = load(file2)["dt"]

bitvec1 = repeat(-.5 .+ (1:size(Icont1)[1]),1,size(Icont1)[2])
lagvec1 = repeat(1:size(Icont1)[2],1,size(Icont1)[1])'

bitvec2 = repeat(-.5 .+ (1:size(Icont2)[1]),1,size(Icont2)[2])
lagvec2 = repeat(1:size(Icont2)[2],1,size(Icont2)[1])'

ylabelsign = [L"sign"]
ylabelexp = [LaTeXString('$'*"e_$i"*'$') for i in 1:8]
ylabelfrac = [LaTeXString('$'*"f_{$i}"*'$') for i in 1:23]

ylabels = cat(ylabelsign,ylabelexp,ylabelfrac,dims=1)

ioff()
fig,(ax1,ax2) = subplots(1,2,figsize=(10,5),sharex=true)
tight_layout(rect=[0.,.02,0.92,0.98])
fig[:subplots_adjust](wspace=0.11,hspace=0.03)
pos2 = ax2[:get_position]()
cax = fig[:add_axes]([pos2[:x1]+0.01,pos2[:y0],0.02,pos2[:y1]-pos2[:y0]])

Q1 = ax1[:pcolormesh](lagvec1,bitvec1,lIcont1,vmin=-3.,vmax=0,cmap="cubehelix_r")
Q2 = ax2[:pcolormesh](lagvec2,bitvec2,lIcont2,vmin=-3.,vmax=0,cmap="cubehelix_r")

ax1[:invert_yaxis]()
ax1[:set_yticks](Array(1:size(Icont1)[1]))
ax1[:set_yticklabels](ylabels)
cb = colorbar(cax=cax,Q2)
cb[:set_label](L"$\log_{10}(I_b)$")

ax2[:invert_yaxis]()
ax2[:set_yticks](Array(1:size(Icont1)[1]))
ax2[:set_yticklabels](ylabels)

tlabels = [1,2,3,5,10,20]
xtiks = [continuous_idx(lags1*dt1/3600/24,t) for t in tlabels]
ax1[:set_xticks](xtiks)
ax1[:set_xticklabels]([latexstring(t) for t in tlabels])
ax1[:set_xlabel]("forecast time [days]")
ax2[:set_xlabel]("forecast time [days]")

ax1[:plot]([0,length(lagvec1)],[1.5,1.5],"grey",lw=1.5)
ax1[:plot]([0,length(lagvec1)],[9.5,9.5],"grey",lw=1.5)

ax2[:plot]([0,length(lagvec2)],[1.5,1.5],"grey",lw=1.5)
ax2[:plot]([0,length(lagvec2)],[9.5,9.5],"grey",lw=1.5)

ax1[:set_ylim](32.5,.5)
ax2[:set_ylim](32.5,.5)
ax1[:set_xlim](1,32)
ax2[:set_xlim](1,32)
ax1[:set_title]("Information content: "*L"u_{i,j}")
ax2[:set_title]("Information content: "*L"u_{i+4,j}")


savefig("figs/juls_inf_cont_uu4.pdf")
close(fig)
