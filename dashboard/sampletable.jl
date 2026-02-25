using YAML

using DataFrames, DataFramesMeta
function metadatafields(sampletable)
    fields = String[]
    for (sample, data) in sampletable["samples"]
        for (entry, value) in data["metadata"]
            push!(fields, entry)
        end
    end
    sort(unique(fields))
end

function selectfirefile(df)
    fire_cram_ind = df.FileClass .== "fire_cram"
    if any(fire_cram_ind)
        return (; selected=fire_cram_ind)
    else
        ### return the fire bam
        if df.Pipeline[1] == "dorado_v13_hg19"
            haplotagged_cram_ind = df.FileClass .== "haplotagged_cram"
            if any(haplotagged_cram_ind)
                return (; selected=haplotagged_cram_ind)
            end
        end
        fire_bam_ind = df.FileClass .== "fire_bam"
        return (; selected=fire_bam_ind)
    end
end

emptymetacol() = Vector{Union{Missing,String}}()

function load_sampletable(sampletable_path::String)


    sampletable = YAML.load_file(sampletable_path)

    @show length(sampletable["samples"])

    sampletable["samples"] = filter(((k, v),) -> get(v, "include", true), sampletable["samples"])


    fields = metadatafields(sampletable)
    samples = sampletable["samples"]
    samplesnames = sort(collect(keys(samples)))



    sampledf = [DataFrame(Sample=samplesnames) DataFrame(mapreduce(f -> [get(samples[s]["metadata"], f, missing) for s in samplesnames], hcat, fields), fields)]
    for n in names(sampledf)
        sampledf[!, n] = identity.(sampledf[!, n])
    end

    mandatory_fields = ["study", "Sample", "cell_type", "mtase", "buffer", "incubation", "num_cells", "adaptive_sampling"]

    sampledf = [sampledf[!, mandatory_fields] sampledf[!, setdiff(names(sampledf), mandatory_fields)]]
    sampledf.adaptive_sampling = coalesce.(sampledf.adaptive_sampling, false)

    filedf = DataFrame(Sample=String[], Pipeline=String[], Genome=String[], FileClass=String[], File=String[])

    for (sample, data) in sampletable["samples"]
        for (pipeline, files) in data["files"]
            genome = occursin("hg38", pipeline) ? "hg38" : "hg19"
            for (fileclass, file) in files
                push!(filedf, (sample, pipeline, genome, fileclass, file))
            end
        end
    end

    alldf = innerjoin(sampledf, filedf, on=:Sample)

    transform!(groupby(alldf, [:Sample, :Pipeline]), selectfirefile)
    alldf.File = joinpath.(sampletable["datadir"], alldf.File)

    alldf
end

# alldf = load_sampletable("test/samples.yaml");