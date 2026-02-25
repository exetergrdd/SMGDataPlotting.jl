"""
    SMGDataPlotting

A package for plotting single molecule genomics data, including fiber-seq and direct RNA-seq from nanopore.
It provides tools for visualizing read pileups, MSPs (Methylation Specific Patches) / FIRE (Fiber-seq inferred regulatory elements) elements, and other genomic data.
"""
module SMGDataPlotting

using SMGReader
using DataFrames, DataFramesMeta, CSV
using CairoMakie

using GenomicFeatures

export pilereads, binvector, firepile, collectfire, assignlevels, loadbed, bedintervals, bedplot, loadgenemodels, filterreads, combinereads, plotgenemodels, fireplot, modplot

include("data.jl")
include("utils.jl")
include("genome.jl")
include("stencillingplots.jl")



end
