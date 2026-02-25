
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

"""
    loadbed(bedfile)

Loads a BED file into a DataFrame. Handles standard BED files and those with headers.
If no header is present, columns are named `:chrom`, `:start`, `:stop`, and optionally `:name`.

# Arguments
- `bedfile`: Path to the BED file.

# Returns
A `DataFrame` containing the BED file data.
"""
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

"""
    bedintervals(bdf)

Converts a DataFrame of BED data (from `loadbed`) into a `GenomicIntervalCollection`.

# Arguments
- `bdf`: DataFrame containing BED data (must have `:chrom`, `:start`, `:stop` columns).

# Returns
A `GenomicIntervalCollection`.
"""
function bedintervals(bdf)
    sort!(bdf, [:chrom, :start, :stop])
    giv = @with bdf GenomicInterval.(:chrom, :start, :stop, '.', :name)
    GenomicIntervalCollection(giv)
end


exonmetadata(geneid, genename, starts, stops) = (; geneid, genename, starts, stops)

"""
    loadgenemodels(exondatafile; filter=:longest)

Loads gene models (exon data) from a file.

# Arguments
- `exondatafile`: Path to the file containing exon data.
- `filter`: Filter strategy (default `:longest`). Currently only supports selecting the longest isoform per gene.

# Returns
A `GenomicIntervalCollection` representing the gene models.
"""
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