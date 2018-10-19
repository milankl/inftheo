using NetCDF
using Printf
using StatsBase
using JLD2
using FileIO

include("infcont.jl")

cd("/home/kloewer/julia/inftheo/")
path = "/local/kloewer/julsdata/"
#path = "/network/aopp/cirrus/pred/kloewer/julsdata/"
runs = [1]

# LOAD DATA
ncu = NetCDF.open(path*"run"*@sprintf("%04d",runs[1])*"/u.nc")
#ncv = NetCDF.open(path*"run"*@sprintf("%04d",runs[1])*"/v.nc")
#nceta = NetCDF.open(path*"run"*@sprintf("%04d",runs[1])*"/eta.nc")

N = size(ncu.vars["u"])[end]
predictor = reshape(ncu.vars["u"][5,50,:],N)
#predictand = reshape(ncu.vars["u"][5,50,:],N)
predictand = predictor

suffix = "uuh"
println(suffix)
dt = ncu.gatts["output_dt"]

# preallocate
lags = cat([0,1,2,3,4,5],Int.(round.(10 .^(0.8:0.06:4.0))),dims=1)
h = 0.02
Icont = infcont(predictor,predictand,h,lags)

save("data/juls/infcont_floats_$suffix.jld2","Icont",Icont,"lags",lags,"dt",dt,"h",h)
