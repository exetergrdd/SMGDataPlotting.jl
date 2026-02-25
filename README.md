# SMGDataPlotting

<!-- [![Build Status](https://github.com/owensnick/SMGDataPlotting.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/owensnick/SMGDataPlotting.jl/actions/workflows/CI.yml?query=branch%3Amain) -->

**SMGDataPlotting.jl** is a Julia package for visualizing single molecule genomics data, such as Fiber-seq and direct RNA-seq from Nanopore. It provides tools to generate high-quality plots for read pileups, chromatin accessibility (FIRE elements), and modified bases.

## Installation

This package is currently in development. You can install it using the Julia package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/owensnick/SMGDataPlotting.jl")
```

## Data Support

SMGDataPlotting works with:
- **BAM/CRAM files**: Support for reading these files is provided by [SMGReader.jl](https://github.com/exetergrdd/SMGReader.jl).
- **BED files**: For genomic intervals and annotations.

## Usage

### Plotting FIRE Elements (Chromatin Accessibility)

The `firepile` function creates a pileup plot of Fiber-seq inferred regulatory elements (FIRE), showing nucleosomes and linkers/MSPs.

```julia
using SMGDataPlotting
using CairoMakie

file = "path/to/your/data.bam"
chrom = "chr1"
loc = 1_000_000:1_005_000

# Create a pileup plot
fig = firepile(file, chrom, loc)
display(fig)
```

### Plotting BED Intervals

You can plot genomic intervals from a BED file using `bedplot`.

```julia
using SMGDataPlotting
using CairoMakie

# Load BED data
bdf = loadbed("path/to/annotations.bed")
intervals = bedintervals(bdf)

# Plot
fig = Figure()
ax = Axis(fig[1, 1])
bedplot("chr1", 1_000_000:1_005_000, ["Region 1"], [intervals]; f=fig, ax=ax)
display(fig)
```

### Collecting Data

You can also collect raw data for custom analysis or plotting:

```julia
# Collect FIRE elements
data = collectfire(file, chrom, loc)
# Access reads and levels
first(data.reads) 

# Collect modified bases and FIRE elements
mod_data = collectmodsfire(file, chrom, loc)
```
