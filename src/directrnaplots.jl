# using SMGReader


const DIRECT_RNA_MODS = (mod_2OmeA, mod_2OmeC, mod_2OmeG, mod_2OmeU, mod_6mA, mod_5mC, mod_inosine, mod_pseU)

const MOD_LABELS = Dict(mod_6mA => "m6A", mod_5mC => "m5C", mod_inosine => "inosine", mod_pseU => "pseU",
    mod_2OmeA => "2Ome-A", mod_2OmeC => "2Ome-C", mod_2OmeG => "2Ome-G", mod_2OmeU => "2Ome-U")


# file = "/Users/ndlo201/projects/directrna/dorado52/12252_RNA/basecall.calledmods.genome.hg38.sort.bam"

# ad = load_directrna_readdata(file)
# length(ad[2])

# chrom, loc = "chr13", 27918000:27928313
# chrom, loc = "chr20", 22578998:22587490
# drd = collectdirectrna(file, chrom, loc)

@inline function incpile!(a, pile, start, off=1)
    i = (a + off) - start + 1
    @inbounds (1 <= i <= length(pile)) && (pile[i] += 1)
    pile
end

@inline function readextentpile!(pile, r, rdata, loc)
    lp = leftposition(r)
    rp = rightposition(rdata)
    for a in lp:rp
        incpile!(a, pile, loc.start)
    end
end


@inline function alignmappile!(pile, rdata, loc)
    for a in rdata.alignmap
        !iszero(a) && incpile!(a, pile, loc.start)
    end
end


function directrnastackplot(chrom, loc, readdata, leveldata; f=Figure(), axv=1, axh=1, ax=Axis(f[axv, axh], xgridvisible=false, ygridvisible=false, xticks=([], [])), exloc=loc)


    readpile = zeros(UInt16, length(exloc))
    modpile = zeros(UInt16, length(exloc), length(DIRECT_RNA_MODS))

    @inbounds for rd in readdata

        for ab in rd.alignblocks
            lp = max(1, ab.start + 1 - first(exloc) + 1)
            rp = min(length(exloc), ab.stop - first(exloc) + 1)
            readpile[lp:rp] .+= 1
        end

        for (k, m) in enumerate(DIRECT_RNA_MODS)
            if haskey(rd.moddata, m)
                for pos in rd.moddata[m]
                    if first(exloc) <= pos <= last(exloc)
                        modpile[pos-first(exloc)+1, k] += 1
                    end
                end
            end
        end
    end

    cc = Makie.wong_colors()

    band!(ax, exloc, 0, readpile, color=:lightgrey)
    for (k, m) in enumerate(DIRECT_RNA_MODS)
        # barplot!(ax, loc, modpile[:, k], color=k > length(cc) ? :black : cc[k], label=string(m), width=20)
        color = k > length(cc) ? :black : cc[k]
        band!(ax, exloc, 0, modpile[:, k], color=color, label=MOD_LABELS[m], strokewidth=1, strokecolor=color)
    end

    hidespines!(ax, :l, :t, :r, :b)
    ymax = 1.1 * maximum(readpile)
    ylims!(ax, 0.0, ymax)
    axislegend(ax, framevisible=false, position=:rc)
    xlims!(ax, first(loc), last(loc))
    f

end



function directrna_blockplot(chrom, loc, readdata, leveldata; f=Figure(), axv=1, axh=1, ax=Axis(f[axv, axh], xgridvisible=false, ygridvisible=false, xticks=([], [])), ymax=1.1 * maximum(leveldata), exloc=loc)


    readrrect = mapreduce((l, read) -> [Rect(ab.start, l + (1 - 0.6) / 2, ab.stop - ab.start, 0.6) for ab in read.alignblocks], vcat, leveldata, readdata)
    !isempty(readrrect) && poly!(ax, readrrect, color=:lightgrey, label="Read")

    cc = Makie.wong_colors()
    w = 1
    for (k, m) in enumerate(DIRECT_RNA_MODS)
        modrect = mapreduce((l, read) -> [Rect(pos - w / 2, l + (1 - 0.9) / 2, w, 0.9) for pos in get(read.moddata, m, Int[])], vcat, leveldata, readdata)


        color = (k > length(cc)) ? :black : cc[k]
        !isempty(modrect) && poly!(ax, modrect, color=color, label=MOD_LABELS[m], strokecolor=color, strokewidth=0.5)
    end

    hidespines!(ax, :l, :t, :r, :b)
    ylims!(ax, 0.0, ymax)
    axislegend(ax, framevisible=false, position=:rc)
    xlims!(ax, first(loc), last(loc))
    f

end

# drd_block = collectdirectrnaablock(file, chrom, loc)

# drd_block[1].moddata
# levels = assignlevels(drd_block);


# @time begin
#     f = Figure()
#     directrna_blockplot(chrom, loc, ad[2], drd_block, levels, f=f, axv=1, axh=1)
#     f
# end
