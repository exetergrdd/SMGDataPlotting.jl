
function bedhasheader(bedfile)
    io = open(bedfile)
    hasheader = false
    while !eof(io)
        line = readline(io)

        if startswith(line, "chrom\tstart\tstop") || startswith(line, "#chrom\tstart\tstop")
            hasheader = true
            break
        elseif startswith(line, "#")
            continue
        else
            break
        end

    end
    close(io)
    hasheader
end

function loadbed(bedfile)
    if bedhasheader(bedfile)
        bdf = CSV.read(bedfile, DataFrame)
        rename!(bdf, replace.(names(bdf), "#chrom" => "chrom"))
    else
        bdf = CSV.read(bedfile, DataFrame, header=false)
        rename!(bdf, 1 => :chrom, 2 => :start, 3 => :stop)
        if size(bdf, 2) >= 4
            rename!(bdf, 4 => :name)
        else
            bdf[!, :name] = string.(replace(basename(bedfile), ".bed" => ""), "_", 1:size(bdf, 1))
        end
    end

    bdf
end

function bedintervals(bdf)
    sort!(bdf, [:chrom, :start, :stop])
    giv = @with bdf GenomicInterval.(:chrom, :start, :stop, '.', :name)
    GenomicIntervalCollection(giv)
end


exonmetadata(geneid, genename, starts, stops) = (; geneid, genename, starts, stops)

function loadgenemodels(exondatafile; filter=:longest)

    exondata = CSV.read(exondatafile, DataFrame)

    if filter == :longest
        exondata = combine(groupby(exondata, [:GeneID, :GeneName]), df -> df[argmax(df.stop .- df.start), :])
    end

    exondata.starts = [parse.(Int, split(s, r"[\[, \]]", keepempty=false)) for s in exondata.starts]
    exondata.stops = [parse.(Int, split(s, r"[\[, \]]", keepempty=false)) for s in exondata.stops]

    givs = @with exondata GenomicInterval.(:chrom, :start, :stop, first.(:strand), exonmetadata.(:GeneID, :GeneName, :starts, :stops))
    GenomicIntervalCollection(givs, true)
end