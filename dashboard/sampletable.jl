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
    fire_cram_ind = df.FileClass .∈ Ref(Set(["fire_cram", "haplotagged_cram", "epi2me_haplotagged_cram", "epi2me_haplotagged_bam", "epi2me_sv_haplotagged_cram", "epi2me_sv_haplotagged_bam"]))
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



function load_directrna_sampletable(sampletable_path::String)

    sampletable = YAML.load_file(sampletable_path)
    @show length(sampletable["samples"])

    sampletable["samples"] = filter(((k, v),) -> get(v, "include", true), sampletable["samples"])
    @show length(sampletable["samples"])

    fields = metadatafields(sampletable)
    samples = sampletable["samples"]
    samplesnames = sort(collect(keys(samples)))

    @show fields
    @show samplesnames
    samples

    sampledf = [DataFrame(Sample=samplesnames) DataFrame(mapreduce(f -> [get(samples[s]["metadata"], f, missing) for s in samplesnames], hcat, fields), fields)]
    for n in names(sampledf)
        sampledf[!, n] = identity.(replace(sampledf[!, n], "N/A" => missing))
    end
    sampledf.pcw = [ismissing(pcw) ? missing : parse(Int, pcw) for pcw in sampledf.pcw]
    sampledf.rin = [ismissing(rin) ? missing : parse(Float64, rin) for rin in sampledf.rin]
    sampledf.tapestation = [ismissing(tapestation) ? missing : parse(Int, tapestation) for tapestation in sampledf.tapestation]
    sampledf.direct_run .= "Y"
    sampledf.Study .= "FetalPancreas"

    mandatory_fields = ["Study", "Sample", "pcw", "sex", "rin", "tapestation"]

    sampledf = [sampledf[!, mandatory_fields] sampledf[!, setdiff(names(sampledf), mandatory_fields)]][!, Not(:direct_run)]

    filedf = DataFrame(Sample=Int[], Pipeline=String[], Genome=String[], FileClass=String[], File=String[])

    for (sample, data) in sampletable["samples"]
        for (pipeline, files) in data["files"]
            genome = data["genome_build"]
            for (fileclass, file) in files
                push!(filedf, (sample, pipeline, genome, fileclass, file))
            end
        end
    end
    filedf = @subset(filedf, :FileClass .∈ Ref(Set(["genome_hg38_sort_bam"])))

    alldf = innerjoin(sampledf, filedf, on=:Sample)
    alldf.File = joinpath.(sampletable["datadir"], alldf.File)

    # push!(alldf, ("FetalPancreas", 2, 1, "M", missing, missing, missing, "pipe", "hg38", "test", "/Users/ndlo201/projects/directrna/dorado52/12252_RNA/basecall.calledmods.genome.hg38.sort.bam"))
    sort!(alldf, [:pcw, :sex])

    alldf
end

# load_directrna_sampletable("test/directrna_samples.yaml")