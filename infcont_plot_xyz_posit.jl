using PyPlot
using JLD
using PyCall
@pyimport cmocean.cm as cmaps

cd("/home/kloewer/julia/inftheo/")

# load one file to get size
Ixx = load("data/infcont_posits_xx.jld")["Icont"]
n,m = size(Ixx)

# preallocate
Iall = zeros(3,3,n,m)

# load all files
for (ip1,predictor) in enumerate(["x","y","z"])
    for (ip2,predictand) in enumerate(["x","y","z"])
        suffix = predictor*predictand
        I = load("data/infcont_posits_$suffix.jld")["Icont"]
        I[I.<=0] = NaN
        Iall[ip1,ip2,:,:] = log10.(I)     # plot information content logarithmically
    end
end

lags = cat(1,[0,1,2,3,4,5],Int.(round.(10.^(-1.5:0.1:2)./0.005)))

bitvec = repmat(-.5 + (1:n),1,m)
lagvec = repmat(1:m,1,n)'

ylabelsign = [L"sign"]
ylabelposit = [LaTeXString('$'*"b_{$i}"*'$') for i in 2:32]
ylabel2 = cat(1,ylabelsign,ylabelposit)

fig,axs = subplots(3,3,figsize=(10,10),sharex=true,sharey=true)
Q = []
for i = 1:3
    for j = 1:3
        Q = axs[i,j][:pcolormesh](lagvec,bitvec,Iall[i,j,:,:],vmin=-4,vmax=0,cmap="cubehelix_r")
    end
end

tight_layout(rect=[0.02,.09,1,0.98])
fig[:subplots_adjust](wspace=0.03,hspace=0.03)

# colorbar
pos = axs[3,1][:get_position]()
pos2 = axs[3,3][:get_position]()
cax = fig[:add_axes]([pos[:x0],0.05,pos2[:x1]-pos[:x0],0.02])
cb = colorbar(cax=cax,Q,orientation="horizontal")
cb[:set_label](L"Information content $\log_{10}(I_b)$")


axs[1,1][:invert_yaxis]()
axs[1,1][:set_yticks](Array(1:n))
axs[1,1][:set_yticklabels](ylabel2,fontsize=7)
axs[2,1][:set_yticklabels](ylabel2,fontsize=7)
axs[3,1][:set_yticklabels](ylabel2,fontsize=7)

# 0.005 is the timestep
xtiks = [find(lags*0.005 .== x) for x in [0.01,0.1,1.,10.,100.]]
axs[3,1][:set_xticks](xtiks)
axs[3,1][:set_xticklabels]([L"10^{-2}",L"10^{-1}",L"10^{0}",L"10^{1}",L"10^{2}"])
axs[3,1][:set_xlabel]("forecast time (mtu)")
axs[3,2][:set_xlabel]("forecast time (mtu)")
axs[3,3][:set_xlabel]("forecast time (mtu)")

#axs[3,1][:plot]([0,36],[1.5,1.5],"w",lw=1.5)
#axs[3,1][:plot]([0,36],[9.5,9.5],"w",lw=1.5)

#ax2[:plot]([0,36],[1.5,1.5],"w",lw=1.5)
#ax2[:plot]([0,36],[9.5,9.5],"k",lw=2)

#ax1[:plot]([7,7],[0,42],"w",lw=1)
#ax2[:plot]([7,7],[0,42],"w",lw=1)

axs[1,1][:set_ylim](32.5,.5)
axs[1,1][:set_xlim](1,43)

subtitl = split("abcdefghi","")
isubtit = 0
for i = 1:3
    for j = 1:3
        axs[i,j][:plot]([0,42],[1.5,1.5],"w",lw=1.5)
        isubtit += 1
        axs[i,j][:text](3,31,subtitl[isubtit],fontweight="bold")
    end
end

axs[1,1][:set_title]("Predictand x")
axs[1,2][:set_title]("Predictand y")
axs[1,3][:set_title]("Predictand z")

axs[1,1][:set_ylabel]("Predictor x",fontsize=12)
axs[2,1][:set_ylabel]("Predictor y",fontsize=12)
axs[3,1][:set_ylabel]("Predictor z",fontsize=12)

savefig("inf_cont_xyz_posits.pdf")
close(fig)
