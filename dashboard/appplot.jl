
mutable struct DisplayTrack{T,C}
    sampleindex::Int
    tracktype::T
    row::Observable{Int}
    file::String
    bystrand::Bool
    byhap::Bool
    plotconfig::C
end
DisplayTrack(sampleindex, tracktype, row, file) = DisplayTrack(sampleindex, tracktype, row, file, false, false, nothing)

function constructaxis(figure, row, col; datatype=:data, config, common_ax_attrs=(xgridvisible=false, ygridvisible=false), internaltitle="", kwargs...)
    ax = Axis(figure[row, col]; common_ax_attrs..., kwargs...)
    rowsize!(figure.layout, row, config.axisheight[datatype])
    ax.yzoomlock = true
    ax.ypanlock = true
    if !isempty(internaltitle)
        text!(ax, 0, 1, text=internaltitle, font=:bold, align=(:left, :top), offset=(4, -2), space=:relative, fontsize=14)
    end
    return ax
end

struct DataCacheKey
    sampleindex::Int
    chrom::String
    loc::UnitRange{Int}
end

struct LevelKey
    hap::Bool
    strand::Bool
end

extend(iv, p=2, l=length(iv), ex=Int(round(p * l))) = (max(1, iv.start .- ex)):(iv.stop.+ex)


@inline function dataplot(chrom, loc, frs, levels, tracktype, hap, strand, ax)

    @show chrom, loc, tracktype, hap, strand
    if tracktype == :FIRE

        #### here implement fireplot with levels
        fireplot(chrom, loc, frs, levels; hap=hap, strand=strand, ax=ax)
    else

        modplot(chrom, loc, frs, levels; mod=tracktype, hap=hap, strand=strand, ax=ax)
    end
end

function browserplot(chrom, loc, displaytracks, sampletable; genemodels=nothing, datacache=Dict{DataCacheKey,Vector{SMGDataPlotting.ModFireRead}}())

    exloc = extend(loc, 2) # pre-extend for fetching data
    @show chrom, loc, exloc, typeof(genemodels)

    axes = Dict{Int,Axis}()
    leveldict = Dict{DataCacheKey,Dict{LevelKey,Vector{Int}}}()
    k = 1
    ## collect data
    for dt in displaytracks
        key = DataCacheKey(dt.sampleindex, chrom, loc)
        if !haskey(datacache, key)
            println("Collecting data for $(sampletable[dt.sampleindex, :Sample]): $key")
            datacache[key] = SMGDataPlotting.collectmodsfire(dt.file, chrom, exloc, byhap=dt.byhap, bystrand=dt.bystrand)
            leveldict[key] = Dict{LevelKey,Vector{Int}}()
        end
        levelkey = LevelKey(dt.byhap, dt.bystrand)
        if !haskey(leveldict[key], levelkey)
            leveldict[key][levelkey] = assignlevels(datacache[key], haplotype=dt.byhap, strand=dt.bystrand)
        end
    end
    config = (size=(1000, 1000), axisheight=Dict(:data => Fixed(200), :genemodel => Fixed(50)))
    figure = Figure(size=config.size)

    #### assess plot rows where each display track has a plot row
    ### each row corresponds to an axis apart from when a display track splits by multiple haplotypes or strands
    ### in this case multiple axes are added, and subsequent rows appear after

    rowdict = Dict{Int,Vector{DisplayTrack}}()
    for dt in displaytracks
        haskey(rowdict, dt.row[]) ? push!(rowdict[dt.row[]], dt) : rowdict[dt.row[]] = [dt]
    end
    # rowhapstrand = Dict{Int,Tuple{Bool,Bool}}()
    rowhap = Dict{Int,Vector{UInt8}}()
    rowstrand = Dict{Int,Vector{Bool}}()
    for (r, dts) in rowdict
        rowhap[r] = UInt8[]
        rowstrand[r] = Bool[]

        for dt in dts

            key = DataCacheKey(dt.sampleindex, chrom, loc)
            frs = datacache[key]
            strands = sort(unique([fr.strand for fr in frs]))
            haplotypes = sort(unique([fr.haplotype for fr in frs]))

            if dt.byhap
                rowhap[r] = union(rowhap[r], haplotypes)
            end
            if dt.bystrand
                rowstrand[r] = union(rowstrand[r], strands)
            end
        end
    end
    rowkeys = sort(collect(keys(rowdict)))
    ### assign axes ranges
    rowaxranges = Dict{Int,UnitRange{Int}}()
    for r in rowkeys
        if !isempty(rowhap[r]) && !isempty(rowstrand[r])
            numaxes = length(rowhap[r]) * length(rowstrand[r])
        elseif !isempty(rowhap[r]) || !isempty(rowstrand[r])
            numaxes = length(rowhap[r]) + length(rowstrand[r])
        else
            numaxes = 1
        end
        rowaxranges[r] = k:(k+numaxes-1)
        k += numaxes
    end


    totalheight = 0
    for dt in displaytracks
        key = DataCacheKey(dt.sampleindex, chrom, loc)
        frs = datacache[key]
        # strands = sort(unique([fr.strand for fr in frs]))
        # haplotypes = sort(unique([fr.haplotype for fr in frs]))
        sample = sampletable[dt.sampleindex, :Sample]
        axesrange = rowaxranges[dt.row[]]
        axc = 0
        if !isempty(rowhap[dt.row[]]) && !isempty(rowstrand[dt.row[]])
            ### add four axes
            for s in rowstrand[dt.row[]], h in rowhap[dt.row[]]
                levels = leveldict[key][LevelKey(true, true)]
                row = axesrange[1] + axc
                ax = constructaxis(figure, row, 1, datatype=:data, config=config, internaltitle=string(sample, "  Tracktype: ", dt.tracktype, " Haplotype: ", h, " Strand: ", s))
                totalheight += config.axisheight[:data].x
                dataplot(chrom, loc, frs, levels, dt.tracktype, h, s, ax)
                axes[row] = ax
                axc += 1
            end


        elseif !isempty(rowstrand[dt.row[]])
            ### add two axes
            for s in rowstrand[dt.row[]]
                levels = leveldict[key][LevelKey(false, true)]
                row = axesrange[1] + axc
                ax = constructaxis(figure, row, 1, datatype=:data, config=config, internaltitle=string(sample, " Tracktype: ", dt.tracktype, " Strand: ", s))
                totalheight += config.axisheight[:data].x
                dataplot(chrom, loc, frs, levels, dt.tracktype, nothing, s, ax)
                axes[row] = ax
                axc += 1
            end

        elseif !isempty(rowhap[dt.row[]])
            ### add two axes
            for h in rowhap[dt.row[]]
                levels = leveldict[key][LevelKey(true, false)]
                row = axesrange[1] + axc
                ax = constructaxis(figure, row, 1, datatype=:data, config=config, internaltitle=string(sample, "  Tracktype: ", dt.tracktype, " Haplotype: ", h))
                totalheight += config.axisheight[:data].x
                dataplot(chrom, loc, frs, levels, dt.tracktype, h, nothing, ax)

                axes[row] = ax
                axc += 1
            end

        else

            row = axesrange[1]
            levels = leveldict[key][LevelKey(false, false)]
            ax = constructaxis(figure, row, 1, datatype=:data, config=config, internaltitle=string(sample, "  Tracktype: ", dt.tracktype))
            totalheight += config.axisheight[:data].x
            dataplot(chrom, loc, frs, levels, dt.tracktype, nothing, nothing, ax)
            axes[row] = ax
        end
    end
    linkyaxes!(values(axes)...)

    if !isnothing(genemodels)
        ax = constructaxis(figure, length(axes) + 1, 1, datatype=:genemodel, config=config, internaltitle="", yreversed=true)
        f, maxlevel = plotgenemodels(chrom, exloc, genemodels; f=figure, ax=ax, vert=true)
        rh = Fixed(maxlevel * 15 + 35)

        rowsize!(figure.layout, length(axes) + 1, rh)
        totalheight += rh.x
        axes[length(axes)+1] = ax
    end

    for (k, ax) in axes
        vlines!(ax, [first(exloc), last(exloc)], color=:black, linestyle=:dash)
        xlims!(ax, first(loc), last(loc))
        if k == 1
            ax.title = string(chrom, ":", first(loc), "-", last(loc))
        end
        if k < length(axes)
            ax.xticks = ([], [])
            rowgap!(figure.layout, k, 8)
        end

    end
    linkxaxes!(values(axes)...)
    # rows = sort(collect(keys(axes)))


    new_figure_height = totalheight + (length(axes) - 1) * 8 + 100

    @show new_figure_height
    resize!(figure.scene, (config.size[1], Int(new_figure_height)))
    figure
end
