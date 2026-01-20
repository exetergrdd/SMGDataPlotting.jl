


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


# function plot