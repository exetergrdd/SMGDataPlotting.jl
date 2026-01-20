module SMGDataPlotting

using SMGReader
using DataFrames, DataFramesMeta, CSV
using CairoMakie

using GenomicFeatures

export  pilereads, binvector, firepile, collectfire, assignlevels, loadbed, bedintervals, bedplot, loadgenemodels

include("data.jl")
include("utils.jl")
include("genome.jl")
include("stencillingplots.jl")



end
