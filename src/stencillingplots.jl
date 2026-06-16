


# plotconfig = (haplo=true, size=(1200, 1000), 
#         tracks=([]),
#         trackdx=10,
#         trackxticks=false,
#         trackheight = Dict(:data => Auto(2.5), :summary => Auto(:0.5), :track => Auto(1.0)),
#         trackmultipler=1.0,
#         # trackheight = Dict(:data => Fixed(100), :summary => Fixed(25), :track => Fixed(30)),
#         filter_unknownhap=false,
#         # config = (()),
#         config=( (data=mod_5mC, plotsum=true),
#                 (data=mod_6mA, plotsum=true),
#                 (data=:fire, plotsum=true)),
#         methconfig=(summary_smooth=100, summary_dx=10, markersize=5, markeralpha=1.0, barsum_width=50, marker=:rect))


# function recordaxes(config)
#     axes = Dict{Symbol, Dict{UInt8, Axis}}()
#     for conf in config

#     end
#     Dict(c => Dict{UInt8, }
# end

function plotfibers(chrom, loc, recorddata, reader; plotconfig, fig=Figure(size=plotconfig.size))

    axisdata = Vector{Tuple{Symbol,Axis}}()

    ### plot top tracks
    if !isempty(plotconfig.tracks)
        ###
    end
    ### plot reads
    ###     1. mods
    ###     2. fire
    ### plot genes etc. for bottom elements

    haplofun = plotconfig.haplo ? (r, rdata) -> haplotype(r, rdata, 0x00) : (r, rdata) -> 0x00

    levels = IntervalLevels(haplofun)


    recordaxes = Dict(c.data => Dict{UInt8,Axis}() for c in plotconfig.config)
    if !isempty(recordaxes)
        for record in eachintersection(reader, chrom, loc)
            processread!(record, recorddata)
            hap = haplofun(record, recorddata)
            l = getlevel(record, recorddata, levels)
            for conf in plotconfig.config
                haskey(recordaxes[conf.data], hap) || (recordaxes[conf][hap] = setuprecordaxes(conf, fig))

                ### now plot fire etc. etc.
            end


        end
    end

end

function setuprecordaxes(conf, fig)
    Axis(fig, yreversed=true, xgridvisible=false, ygridvisible=false)
end

@inline function plotfire!(ax, record, recorddata, level)
    firecoordconvert = firegenome(f.r, f.rdata)
    nucs = firecoordconvert(firenucs(record, recorddata))
    msps = firecoordconvert(firemsps(record, recorddata))



    linkerrects = Rect[]
    firerects = Rect[]
    firecols = Float64[]
    for (start, stop, q) in msps
        r = Rect(s)

        if q < 0.9 * 255
            push!(linkerrects, Rect(start, l + (1 - 0.9) / 2, stop - start, 0.9))
        else

        end

    end
    poly!(axd, [Rect(start, l + (1 - 0.9) / 2, stop - start, 0.9) for (start, stop, q) in msps if q < 0.9 * 255], color=:dodgerblue, label="Linker")

    ### nucleosomes
    # for (start, stop) in 

    # end
end




"""
    firepile(file, chrom, loc; f = Figure(), ax=Axis(f[1, 1]))

Creates a pileup plot of FIRE elements (nucleosomes and linkers/MSPs) for a given genomic region.

# Arguments
- `file`: Path to the BAM/CRAM file.
- `chrom`: Chromosome name.
- `loc`: UnitRange representing the genomic locus.
- `f`: A Makie `Figure` object (optional, defaults to a new figure).
- `ax`: A Makie `Axis` object (optional, defaults to `f[1, 1]`).

# Returns
The Makie `Figure` object containing the plot.
"""
function firepile(file, chrom, loc; f=Figure(), ax=Axis(f[1, 1]))
    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapModFire())

    levels = IntervalLevels()


    for record in eachintersection(reader, chrom, loc)
        validflag(record) || continue
        processread!(record, recorddata) || continue
        firegenomecoords = firegenome(record, recorddata)
        l = getlevel(record, recorddata, levels)

        nucs = firegenomecoords(firenucs(record, recorddata))
        poly!(ax, [Rect(start, l + (1 - 0.6) / 2, stop - start, 0.6) for (start, stop) in nucs], color=:lightgrey, label="Nucleosome")

        msps = firegenomecoords(firemsps(record, recorddata))
        poly!(ax, [Rect(start, l + (1 - 0.9) / 2, stop - start, 0.9) for (start, stop, q) in msps if q < 0.9 * 255], color=:dodgerblue, label="Linker")
        poly!(ax, [Rect(start, l + (1 - 0.9) / 2, stop - start, 0.9) for (start, stop, q) in msps if q >= 0.9 * 255], color=:red, label="FIRE")
    end

    close(reader)

    f
end



# fireplot(chrom, loc, firereads; kwargs...) = fireplot(chrom, loc, firereads.reads, firereads.levels; kwargs...)


function fireplot(chrom, loc, readdata, leveldata; hap=nothing, strand=nothing, f=Figure(), axv=1, axh=1, ax=Axis(f[axv, axh], xgridvisible=false, ygridvisible=false, xticks=([], [])), ymax=1.1 * maximum(leveldata), showlegend=true, expand_fire = 1, rasterize=false)


    if !isnothing(hap) && !isnothing(strand)
        ind = [(r.haplotype == hap) && (r.strand == strand) for r in readdata]
        reads = @view readdata[ind]
        levels = @view leveldata[ind]
    elseif !isnothing(hap)
        hapind = [r.haplotype == hap for r in readdata]
        reads = @view readdata[hapind]
        levels = @view leveldata[hapind]
    elseif !isnothing(strand)
        strandind = [r.strand == strand for r in readdata]
        reads = @view readdata[strandind]
        levels = @view leveldata[strandind]
    else
        reads = readdata
        levels = leveldata
    end


    lr = @time mapreduce((l, fr) -> [Rect(start, l + (1 - 0.9) / 2, stop - start, 0.9) for (start, stop, q) in fr.msps if q < 0.9 * 255], vcat, levels, reads)
    mr = @time mapreduce((l, fr) -> [Rect(start - (stop - start) * expand_fire, l + (1 - 0.9) / 2, (stop - start) * (1 + expand_fire), 0.9) for (start, stop, q) in fr.msps if q >= 0.9 * 255], vcat, levels, reads)


    # @show length(nr)
    # @show length(lr)
    # @show length(mr)

    nr = @time mapreduce((l, fr) -> [Rect(start, l + (1 - 0.6) / 2, stop - start, 0.6) for (start, stop) in fr.nucs], vcat, levels, reads)
    !isempty(nr) && poly!(ax, nr, color=:lightgrey, label="Nucleosome", rasterize=rasterize)
    # @show "one nuc"
    # nr = [Rect(fr.nucs[1][1], l + (1 - 0.6) / 2, fr.nucs[end][2] - fr.nucs[1][1], 0.6) for (fr, l) in zip(reads, levels) if !isempty(fr.nucs)]
    # !isempty(nr) && poly!(ax, nr, color=:lightgrey, label="Nucleosome")

    !isempty(lr) && poly!(ax, lr, color=:dodgerblue, label="Linker", rasterize=rasterize)
    !isempty(mr) && poly!(ax, mr, color=:red, label="FIRE", rasterize=rasterize)


    hidespines!(ax, :l, :t, :r, :b)
    ylims!(ax, 0.0, ymax)
    showlegend && axislegend(ax, framevisible=false, position=:rc)
    xlims!(ax, first(loc), last(loc))
    f
end

function modplot(chrom, loc, readdata, leveldata; mod=mod_6mA, hap=nothing, showline=true, validpos=(Set{Int}(), Set{Int}()), strand=nothing, f=Figure(), axv=1, axh=1, ax=Axis(f[axv, axh], xgridvisible=false, ygridvisible=false, xticks=([], [])), ymax=1.1 * maximum(leveldata), markersize=4, markeralpha=1.0)

    if !isnothing(hap) && !isnothing(strand)
        ind = [(r.haplotype == hap) && (r.strand == strand) for r in readdata]
        reads = @view readdata[ind]
        levels = @view leveldata[ind]
    elseif !isnothing(hap)
        hapind = [r.haplotype == hap for r in readdata]
        reads = @view readdata[hapind]
        levels = @view leveldata[hapind]
    elseif !isnothing(strand)
        strandind = [r.strand == strand for r in readdata]
        reads = @view readdata[strandind]
        levels = @view leveldata[strandind]
    else
        reads = readdata
        levels = leveldata
    end

    if Symbol(mod) == :mod_6mA

        modfield = :mods_6mA
        c = :steelblue
    elseif Symbol(mod) == :mod_5mC
        modfield = :mods_5mC
        c = :green
    elseif Symbol(mod) == :mod_5hmC
        modfield = :mods_5hmC
        c = :orange
    end

    for (fr, l) in zip(reads, levels)
        showline && lines!(ax, [fr.lp, fr.rp], [l, l], color=c)
        data = getproperty(fr, modfield)

        if !isempty(validpos[1])

            if fr.strand
                data = filter(d -> d ∈ validpos[1], data)
            else
                data = filter(d -> d ∈ validpos[2], data)
            end
        end
        scatter!(ax, data, fill(l, length(data)), marker=:rect, color=c, markersize=markersize, alpha=markeralpha)
    end


    hidespines!(ax, :l, :t, :r, :b)
    ylims!(ax, 0, ymax)
    # axislegend(ax, framevisible=false, position=:rc)
    xlims!(ax, first(loc), last(loc))
    f
end



function modallplot(chrom, loc, readdata, leveldata; mod=mod_6mA, hap=nothing, strand=nothing, unmod=false, plotline=true, f=Figure(), axv=1, axh=1, ax=Axis(f[axv, axh], xgridvisible=false, ygridvisible=false, xticks=([], [])), ymax=1.1 * maximum(leveldata), markersize=4, markeralpha=1.0, color=nothing)

    if !isnothing(hap) && !isnothing(strand)
        ind = [(r.haplotype == hap) && (r.strand == strand) for r in readdata]
        reads = @view readdata[ind]
        levels = @view leveldata[ind]
    elseif !isnothing(hap)
        hapind = [r.haplotype == hap for r in readdata]
        reads = @view readdata[hapind]
        levels = @view leveldata[hapind]
    elseif !isnothing(strand)
        strandind = [r.strand == strand for r in readdata]
        reads = @view readdata[strandind]
        levels = @view leveldata[strandind]
    else
        reads = readdata
        levels = leveldata
    end

    if Symbol(mod) == :mod_6mA

        modfield = ifelse(unmod, :un_mods_6mA, :mods_6mA)
        c = :steelblue
    elseif Symbol(mod) == :mod_5mC
        # modfield = :mods_5mC
        modfield = ifelse(unmod, :un_mods_5mC, :mods_5mC)
        c = :green
    elseif Symbol(mod) == :mod_5hmC
        # modfield = :mods_5hmC
        modfield = ifelse(unmod, :un_mods_5hmC, :mods_5hmC)
        c = :orange
    end

    if color !== nothing
        c = color
    end

    for (fr, l) in zip(reads, levels)
        plotline && lines!(ax, [fr.lp, fr.rp], [l, l], color=c)
        data = getproperty(fr, modfield)
        if unmod
            scatter!(ax, data, fill(l, length(data)), marker=:rect, color=:white, strokecolor=c, strokewidth=1, markersize=markersize * 0.75, alpha=markeralpha)
        else
            scatter!(ax, data, fill(l, length(data)), marker=:rect, color=c, markersize=markersize, alpha=markeralpha)
        end
    end


    hidespines!(ax, :l, :t, :r, :b)
    ylims!(ax, 0, ymax)
    # axislegend(ax, framevisible=false, position=:rc)
    xlims!(ax, first(loc), last(loc))
    f
end



"""
    bedplot(chrom, loc, labels, ivs; f=Figure(), axh=1, axv=1, ax=...)

Plots BED intervals as polygons.

# Arguments
- `chrom`: Chromosome name.
- `loc`: UnitRange representing the genomic locus.
- `labels`: Vector of labels for the intervals.
- `ivs`: Vector of intervals to plot (e.g., from `readintervals`).
- `f`: Makie `Figure` object.
- `ax`: Makie `Axis` object.

# Returns
The Makie `Figure` object.
"""
function bedplot(chrom, loc, labels, ivs; f=Figure(), axh=1, axv=1, ax=Axis(f[axv, axh], xgridvisible=false, ygridvisible=false, yreversed=true, yticks=([], [])))
    k = 1
    for (l, iv) in zip(labels, ivs)
        poly!(ax, [Rect(ii.first, k + (1 - 0.6) / 2, ii.last - ii.first, 0.6) for ii in eachoverlap(iv, GenomicInterval(chrom, loc))], label=l)
        k += 1
    end
    hidespines!(ax, :t, :l, :r, :b)
    axislegend(ax, framevisible=false)
    ax.xticks = ([], [])
    ylims!(0.5, k + 1)
    f
end


function genefilter(gi)
    data = GenomicFeatures.metadata(gi)
    if occursin(r"^ENS|^AC[0-9]*", data.genename)
        return false
    else
        return true
    end
end


function genetext(chrom, loc, iv)
   
    
    
end

function plotgenemodels(chrom, loc, giv; f=Figure(), axh=1, axv=1, ax=Axis(f[axv, axh], xgridvisible=false, ygridvisible=false), vert=true, genefilt=genefilter, addlevel = Dict{String, Int}())

    k = 1
    off = length(loc) * 0.005

    intervals = [gi for gi in eachoverlap(giv, GenomicInterval(chrom, loc)) if genefilt(gi)]
    

    if isempty(intervals)
        levels = Int[]
        maxlevel = 0
    else
        levels = assignlevels(intervals)
        maxlevel = maximum(levels)
    end

    row_height = 1
    row_spacing = 1.3
    exon_height = 0.6 ### of row_height
    label_space = 0.2 ### of row_height
    for (iv, l) in zip(intervals, levels)
        data = GenomicFeatures.metadata(iv)
        @show data.genename, l, strand(iv)
        l += get(addlevel, data.genename, 0)
        y = l * (row_height + row_spacing)
        lines!(ax, [iv.first, iv.last], [y + row_height / 2, y + row_height / 2], color=:black)
        poly!(ax, [Rect(start, y + row_height * (1 - exon_height) / 2, stop - start, row_height * exon_height) for (start, stop) in zip(data.starts, data.stops)], color=:black)
        
        
        ### mark arrows for direction
        if (strand(iv) == STRAND_POS) || (strand(iv) == STRAND_NEG)
            
            if length(data.stops) == 1
                ap = [(data.starts[1] + data.stops[1])/2]
                color = :white
            else
                ap = zeros(length(data.stops) - 1)
                for i = 1:(length(data.stops) - 1)
                    
                    if strand(iv) == STRAND_POS
                        ap[i] = (data.stops[i] + data.starts[i+1])/2
                    else
                        ap[i] = (data.starts[i] + data.stops[i+1])/2
                    end
                end
                color = :black
            end



            scatter!(ax, ap, fill(y + row_height / 2, length(ap)), marker=strand(iv) == STRAND_POS ? '>' : '<', color=color)
        end
        
        gene_x = div(max(iv.first, first(loc)) + min(iv.last, last(loc)), 2)
        
        text!(ax, gene_x, y + row_height * label_space, text=data.genename, align=(:center, :top), fontsize=12)
        # vert && (k += 1)
    end

    xlims!(ax, first(loc), last(loc))
    hidespines!(ax, :r, :l, :t)
    ax.yticks = ([], [])
    ylims!(ax, 1 + 0.5 - 1, maxlevel * (row_height + row_spacing) + 1)
    f, maxlevel
end