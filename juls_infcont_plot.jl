using PyPlot
using JLD
using PyCall
@pyimport cmocean.cm as cmaps

cd("/home/kloewer/julia/inftheo/")

Icontf = load("data/juls_infcont_floats_uu.jld")["Icont"]
Icontp = load("data/juls_infcont_floats_uu1.jld")["Icont"]

Icontf[Icontf.<=0] = NaN
lIcontf = log10.(Icontf)

Icontp[Icontp.<=0] = NaN
lIcontp = log10.(Icontp)

#lags = Int.(round.(10.^(-1.5:0.1:2)./0.005))
lags = cat(1,[0,1,2,3,4,5],Int.(round.(10.^(0.8:0.06:3))))

bitvec = repmat(-.5 + (1:size(Icontf)[1]),1,size(Icontf)[2])
lagvec = repmat(1:size(Icontf)[2],1,size(Icontf)[1])'

ylabelsign = [L"sign"]
ylabelexp = [LaTeXString('$'*"e_$i"*'$') for i in 1:8]
ylabelfrac = [LaTeXString('$'*"f_{$i}"*'$') for i in 1:23]
ylabelposit = [LaTeXString('$'*"b_{$i}"*'$') for i in 2:32]

ylabels = cat(1,ylabelsign,ylabelexp,ylabelfrac)
ylabel2 = cat(1,ylabelsign,ylabelposit)

ioff()
fig,(ax1,ax2) = subplots(1,2,figsize=(10,5),sharex=true)
tight_layout(rect=[0.,.02,0.92,0.98])
fig[:subplots_adjust](wspace=0.11,hspace=0.03)
pos2 = ax2[:get_position]()
cax = fig[:add_axes]([pos2[:x1]+0.01,pos2[:y0],0.02,pos2[:y1]-pos2[:y0]])

Q1 = ax1[:pcolormesh](lagvec,bitvec,lIcontf,vmin=-2.,vmax=0,cmap="cubehelix_r")
Q2 = ax2[:pcolormesh](lagvec,bitvec,lIcontp,vmin=-2.,vmax=0,cmap="cubehelix_r")

ax1[:invert_yaxis]()
ax1[:set_yticks](Array(1:size(Icontf)[1]))
ax1[:set_yticklabels](ylabels)
cb = colorbar(cax=cax,Q2)
cb[:set_label](L"$\log_{10}(I_b)$")

ax2[:invert_yaxis]()
ax2[:set_yticks](Array(1:size(Icontf)[1]))
ax2[:set_yticklabels](ylabels)

xtiks = [find(lags .== x) for x in [4,8,25.,200.]]
ax1[:set_xticks](xtiks)
ax1[:set_xticklabels]([L"1",L"2",L"6",L"50"])
ax1[:set_xlabel]("forecast time [days]")
ax2[:set_xlabel]("forecast time [days]")

ax1[:plot]([0,length(lagvec)],[1.5,1.5],"grey",lw=1.5)
ax1[:plot]([0,length(lagvec)],[9.5,9.5],"grey",lw=1.5)

ax2[:plot]([0,length(lagvec)],[1.5,1.5],"grey",lw=1.5)
ax2[:plot]([0,length(lagvec)],[9.5,9.5],"grey",lw=1.5)
#ax2[:plot]([0,36],[9.5,9.5],"k",lw=2)

#ax1[:plot]([7,7],[0,42],"w",lw=1)
#ax2[:plot]([7,7],[0,42],"w",lw=1)


ax1[:set_ylim](32.5,.5)
ax2[:set_ylim](32.5,.5)
ax1[:set_xlim](1,32)
ax2[:set_xlim](1,32)
ax1[:set_title]("32bit Floats")
ax2[:set_title]("32bit Floats")


savefig("figs/juls_inf_cont.pdf")
close(fig)
