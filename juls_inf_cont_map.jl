using NetCDF
using Printf
using StatsBase
using JLD2
using FileIO

include("infcont.jl")

cd("/home/kloewer/julia/inftheo/")
#path = "/local/kloewer/julsdata/"
path = "/network/aopp/cirrus/pred/kloewer/julsdata/"
runs = [1]

# LOAD DATA
ncu = NetCDF.open(path*"run"*@sprintf("%04d",runs[1])*"/u2.nc")
#ncv = NetCDF.open(path*"run"*@sprintf("%04d",runs[1])*"/v.nc")
#nceta = NetCDF.open(path*"run"*@sprintf("%04d",runs[1])*"/eta.nc")

N = size(ncu.vars["u"])
suffix = "uu"
dt = ncu.gatts["output_dt"]

# preallocate
Icont = Array{Float64,3}(undef,N[1],N[2],32)

# options
lags = [0]
h = 0.02

for i in 1:2
    for j in 1:2

        # read data
        predictor = reshape(ncu.vars["u"][i,j,1:N[end]],N[end])
        predictand = predictor

        println("Data read.")
        Icont[i,j,:] = infcont(predictor,predictand,h,lags)
        println("$i/$(N[1]), $j/$(N[2])")
    end
end

save("data/juls/infcont_map_$suffix.jld2","Icont",Icont,"lags",lags,"dt",dt,"h",h)
