using PyPlot
using JLD
using PyCall

cd("/home/kloewer/julia/inftheo/")

# load one file to get size
Ixx = load("data/infcont_floats_xx.jld")["Icont"]
n,m = size(Ixx)

# preallocate
Isum = zeros(n,m)

# load all files
for predictor in ["x","y","z"]
    for predictand in ["x","y","z"]
        suffix = predictor*predictand
        Isum += load("data/infcont_floats_$suffix.jld")["Icont"]
    end
end

Isum[Isum.<0] = NaN
Isum = log10.(Isum)

lags = cat(1,[0,1,2,3,4,5],Int.(round.(10.^(-1.5:0.1:2)./0.005)))

bitvec = repmat(-.5 + (1:n),1,m)
lagvec = repmat(1:m,1,n)'

ylabelsign = [L"sign"]
ylabelexp = [LaTeXString('$'*"e_$i"*'$') for i in 1:8]
ylabelfrac = [LaTeXString('$'*"f_{$i}"*'$') for i in 1:23]
ylabelposit = [LaTeXString('$'*"b_{$i}"*'$') for i in 2:32]

ylabels = cat(1,ylabelsign,ylabelexp,ylabelfrac)
ylabel2 = cat(1,ylabelsign,ylabelposit)

fig,axs = subplots(figsize=(6,5))
Q = axs[:pcolormesh](lagvec,bitvec,Isum,vmin=-3,vmax=1,cmap="cubehelix_r")

tight_layout(rect=[0.02,.15,1,0.98])
#fig[:subplots_adjust](wspace=0.03,hspace=0.03)

# colorbar
pos = axs[:get_position]()
pos2 = axs[:get_position]()
cax = fig[:add_axes]([pos[:x0],0.1,pos2[:x1]-pos[:x0],0.02])
cb = colorbar(cax=cax,Q,orientation="horizontal")
cb[:set_label](L"Information content $\log_{10}(I_b)$")


axs[:invert_yaxis]()
axs[:set_yticks](Array(1:size(Icontf)[1]))
axs[:set_yticklabels](ylabels,fontsize=7)

xtiks = [find(lags*0.005 .== x) for x in [0.01,0.1,1.,10.,100.]]
axs[:set_xticks](xtiks)
axs[:set_xticklabels]([L"10^{-2}",L"10^{-1}",L"10^{0}",L"10^{1}",L"10^{2}"])
axs[:set_xlabel]("forecast time (mtu)")

#axs[3,1][:plot]([0,36],[1.5,1.5],"w",lw=1.5)
#axs[3,1][:plot]([0,36],[9.5,9.5],"w",lw=1.5)

#ax2[:plot]([0,36],[1.5,1.5],"w",lw=1.5)
#ax2[:plot]([0,36],[9.5,9.5],"k",lw=2)

#ax1[:plot]([7,7],[0,42],"w",lw=1)
#ax2[:plot]([7,7],[0,42],"w",lw=1)

axs[:set_ylim](32.5,.5)
axs[:set_xlim](1,43)

axs[:plot]([0,42],[1.5,1.5],"w",lw=1.5)
axs[:plot]([0,42],[9.5,9.5],"w",lw=1.5)

axs[:set_title]("Information content (x,y,z)")
savefig("inf_cont_xyz_sum.pdf")
close(fig)
