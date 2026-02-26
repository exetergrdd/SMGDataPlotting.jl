using Bonito
using WGLMakie
using DataFrames
using Glob
using SMGReader
using SMGDataPlotting
using GenomicFeatures

include("sampletable.jl")
include("appplot.jl")


function autodetectbamdata(bampath)
    htsdata = autodetecthtsdata(bampath)
    auxmap = autodetectaux(bampath)
    mods = autodetectmods(bampath)
    hasfire = typeof(auxmap) <: AuxMapModFire
    (; hasfire, mods)
end


function load_annotated_samples(yamlfile="test/samples.yaml")
    alldf = @subset(load_sampletable(yamlfile), :selected)

    alldf.FIRE = falses(size(alldf, 1))
    alldf.Mods = [Modification[] for _ in 1:size(alldf, 1)]

    for i in 1:size(alldf, 1)
        try
            baminfo = autodetectbamdata(alldf.File[i])
            alldf.FIRE[i] = baminfo.hasfire
            alldf.Mods[i] = collect(baminfo.mods)
        catch
            ### silently ignore
        end
    end
    alldf = @subset(alldf, :FIRE .| length.(:Mods) .> 0)
    excludefields = ["selected", "notes"]
    alldf = alldf[!, setdiff(names(alldf), excludefields)]
    alldf
end



function test_load()
    dirs = glob("*/", "/Users/ndlo201/projects/smf/fire/hap/")
    cramfiles = mapreduce(d -> glob("*.cram", joinpath(d, "fire/")), vcat, dirs)
    samples = basename.(dirname.(dirname.(cramfiles)))

    df = DataFrame(Study=String[], Sample=String[], MTase=String[], AdaptiveSample=Bool[], Genome=String[], FIRE=Bool[], Mods=Vector{Modification}[], File=String[])

    for (sample, file) in zip(samples, cramfiles)

        try
            craminfo = autodetectbamdata(file)
            push!(df, ("Hap", sample, "EcoGII", false, "hg19", craminfo.hasfire, collect(craminfo.mods), file))
        catch
            ### silently ignore
        end
    end

    push!(df, ("Hap", "test", "EcoGII", false, "hg38", true, [mod_6mA, mod_5mC], ""))

    return df
end

function loadgenedata(genome="hg19", host=gethostname())
    if host == "penrose"
        if genome == "hg19"
            gene_ivs = loadgenemodels("/penrose/resource/human/hg19/gencode.v45lift37.annotation.exondata.tsv.gz")
        elseif genome == "hg38"
            gene_ivs = loadgenemodels("/penrose/resource/human/hg38/gencode.v43.annotation.exondata.tsv.gz")
        end
    else
        if genome == "hg19"
            gene_ivs = loadgenemodels("/Users/ndlo201/resource/human/hg19/gencode.v45lift37.annotation.exondata.tsv.gz")
        elseif genome == "hg38"
            gene_ivs = loadgenemodels("/Users/ndlo201/resource/human/hg38/gencode.v43.annotation.exondata.tsv.gz")
        end
    end

    (genedict=Dict(GenomicFeatures.metadata(iv).genename => iv for iv in gene_ivs), gene_ivs=gene_ivs)
end

function adaptivedir(host=gethostname())
    if host == "penrose"
        return "/penrose/projects/ont/resource/as_library"
    else
        return "/Users/ndlo201/projects/smf/adaptive_sampling/library"
    end
end

targetgenome(target) = replace(basename(target), "targets." => "", ".bed" => "")

function loadadaptivesamplingregions()
    dirs = glob("*/", adaptivedir())

    samples = basename.(dirname.(dirs))

    targetdict = Dict(s => glob("targets*", d) for (s, d) in zip(samples, dirs))


    genomes = sort(unique(mapreduce(fs -> targetgenome.(fs), vcat, values(targetdict))))

    ast_files = Dict(g => Dict{String,String}() for g in genomes)

    for (s, files) in targetdict
        for f in files
            g = targetgenome(f)
            ast_files[g][s] = f
        end
    end

    astargets = Dict(g =>
        Dict(s =>
            Dict(String(GenomicFeatures.metadata(iv)) =>
                iv for iv in bedintervals(loadbed(f)))
             for (s, f) in d)
                     for (g, d) in ast_files)

    astargets
end


function makedisplaytable(allmetadata, genome)


    metadata = allmetadata[allmetadata.Genome.==genome, :]

    firecheckboxes = [metadata.FIRE[i] ? Bonito.Checkbox(false) : missing for i in 1:size(metadata, 1)]
    mod_6mA_checkboxes = [mod_6mA ∈ metadata.Mods[i] ? Bonito.Checkbox(false) : missing for i in 1:size(metadata, 1)]
    mod_5mC_checkboxes = [mod_5mC ∈ metadata.Mods[i] ? Bonito.Checkbox(false) : missing for i in 1:size(metadata, 1)]
    mod_5hmC_checkboxes = [mod_5hmC ∈ metadata.Mods[i] ? Bonito.Checkbox(false) : missing for i in 1:size(metadata, 1)]



    dtable = deepcopy(metadata[!, Not([:FIRE, :Mods])])
    dtable.FIRE = firecheckboxes
    dtable.mod_6mA = mod_6mA_checkboxes
    dtable.mod_5mC = mod_5mC_checkboxes
    dtable.mod_5hmC = mod_5hmC_checkboxes

    dtable
end





function checkboxlistener(displaytable, label, samples)

    for (k, (cb, file)) in enumerate(zip(displaytable[!, label], displaytable.File))
        ismissing(cb) && continue
        on(cb.value) do val
            if val

                if isempty(samples[])
                    maxrow = 1
                else
                    maxrow = maximum(s -> s.row[], samples[]) + 1
                end
                push!(samples[], DisplayTrack(k, label, Observable(maxrow), file))
            else

                ind = findfirst(s -> (s.sampleindex == k) && (s.tracktype == label), samples[])
                if !isnothing(ind)
                    deleteat!(samples[], ind)
                end

            end
            notify(samples)
        end
    end

end


function parsesearchbox(val, genemodels)
    val = strip(val)
    isempty(val) && return nothing

    val = replace(val, "," => "")

    # Unified coordinate regex
    m = match(r"^([a-zA-Z0-9._]+)[\s:]+(\d+)[\s-]+(\d+)$", val)
    if !isnothing(m)
        chr = m[1]
        start = parse(Int, m[2])
        stop = parse(Int, m[3])
        return chr, start, stop
    end

    # Try gene lookup
    gene = uppercase(val)
    if haskey(genemodels.genedict, gene)
        iv = genemodels.genedict[gene]
        return (
            GenomicFeatures.seqname(iv),
            GenomicFeatures.leftposition(iv),
            GenomicFeatures.rightposition(iv)
        )
    end

    return nothing
end


function as_region_sort_lt(regA, regB)
    mA = match(r"^[A-Za-z0-9]+_\d+", regA)
    mB = match(r"^[A-Za-z0-9]+_\d+", regB)

    if isnothing(mA) && isnothing(mB)
        return regA < regB
    elseif isnothing(mA)
        return true
    elseif isnothing(mB)
        return false
    else
        return regA < regB
    end
end


function smgbrowser(metadata=test_load())


    app = App() do session


        ###### styles

        main_style = Styles(
            CSS(
                "font-family" => "sans-serif",
                "display" => "flex",
                "flex-direction" => "column",
                "gap" => "2px"
            )
        )

        tight_style = Styles(
            CSS(
                "padding" => "4px 8px",
                "font-size" => "0.9em",
                "height" => "auto",
                "width" => "fit-content",
                "min-height" => "20px",
                "line-height" => "1.2",
                "border-radius" => "4px",
                "border" => "1px solid #ccc"
            )
        )

        minimal_style = Styles(
            CSS(
                "padding" => "4px 8px",
                "font-size" => "1em",
                "width" => "fit-content",
                "border" => "none",
                "background" => "#f0f0f0",
                "border-radius" => "4px"
            )
        )

        wide_text_style = Styles(
            CSS(
                "padding" => "4px 8px",
                "font-size" => "1em",
                "width" => "200px",      # ← fixed width
                "border" => "none",
                "background" => "#f0f0f0",
                "border-radius" => "4px"
            )
        )


        #### load samples and targets
        genome = Bonito.Dropdown(["hg19", "hg38"], style=minimal_style)
        samples = Observable(Vector{DisplayTrack}())


        ### load gene models
        defaultgenome = genome.value[]
        genemodelcache = Dict(defaultgenome => loadgenedata(defaultgenome))
        genemodels = Observable(genemodelcache[defaultgenome])

        on(genome.value; update=true) do g
            if !haskey(genemodelcache, g)
                genemodelcache[g] = loadgenedata(g)
            end
            genemodels[] = genemodelcache[g]
        end


        #### adaptive sampling selection and update logic
        astargets = loadadaptivesamplingregions()
        ast_genome = lift(g -> astargets[g], genome.value)
        ast_libraries = lift(g -> sort(collect(keys(g))), ast_genome)

        ast_selected_library = Observable{String}(first(ast_libraries[]))
        # ast_selected_region = Observable{String}("")
        ast_library_dropdown = Bonito.Dropdown([""], style=minimal_style)
        ast_region_dropdown = Bonito.Dropdown([""], style=minimal_style)
        ast_load_button = Bonito.Button("Load", style=tight_style)


        ### on load or update populate library dropdown
        on(ast_genome; update=true) do ag
            libs = sort(collect(keys(ag)))
            ast_library_dropdown.options[] = libs

            if !isempty(libs)
                ast_library_dropdown.value[] = first(libs)
            end
        end
        ### update Observable
        on(ast_library_dropdown.value) do val
            ast_selected_library[] = val
        end
        ### update region options based on selected library
        ast_regions = lift(ast_genome, ast_selected_library) do ag, lib
            sort(collect(keys(ag[lib])), lt=as_region_sort_lt)
        end
        ### on load or update populate region dropdown
        on(ast_regions; update=true) do regs
            ast_region_dropdown.options[] = regs
            ast_region_dropdown.value[] = ""   # no auto-select

            if !isempty(regs)
                ast_region_dropdown.value[] = first(regs)
            end
        end
        ### update Observabl


        #### region navigation
        chrom = Observable("chr1")
        start = Observable(1)
        stop = Observable(1000)

        locationsearch = Bonito.TextField("", style=wide_text_style)

        #### load genes and gene search
        on(locationsearch.value) do val
            parsed = parsesearchbox(val, genemodels[])
            @show parsed
            if !isnothing(parsed)
                chrom.val = parsed[1]
                start.val = parsed[2]
                stop.val = parsed[3]
                notify(chrom)
            else
                @warn "Invalid location search query: $(val)"
            end
        end
        gobutton = Bonito.Button("Go", style=tight_style)
        on(gobutton.value) do val
            parsed = parsesearchbox(locationsearch.value[], genemodels[])
            @show parsed
            if !isnothing(parsed)
                chrom.val = parsed[1]
                start.val = parsed[2]
                stop.val = parsed[3]
                notify(chrom)
            else
                @warn "Invalid location search query: $(locationsearch.value[])"
            end
        end

        on(ast_load_button.value) do val
            reg = ast_region_dropdown.value[]

            gi = ast_genome[][ast_selected_library[]][reg]
            # chrom.val = GenomicFeatures.seqname(gi)
            # start.val = GenomicFeatures.leftposition(gi)
            # stop.val = GenomicFeatures.rightposition(gi)
            locationsearch.value[] = string(GenomicFeatures.seqname(gi), ":", GenomicFeatures.leftposition(gi), "-", GenomicFeatures.rightposition(gi))
        end


        ###### Sample logic

        samplecheckboxtable = lift(genome.value) do g
            empty!(samples[])
            notify(samples)
            dtable = makedisplaytable(metadata, g)
            checkboxlistener(dtable, :FIRE, samples)
            checkboxlistener(dtable, :mod_6mA, samples)
            checkboxlistener(dtable, :mod_5mC, samples)
            checkboxlistener(dtable, :mod_5hmC, samples)
            dtable
        end
        displaytable = lift(d -> Bonito.Table(d[!, Not(:File)]), samplecheckboxtable)


        selectedsamples = lift(samples) do samples

            if isempty(samples)
                return DataFrame()
            end

            inds = [s.sampleindex for s in samples]
            rowtextfield = [TextField(string(s.row[]), style=minimal_style) for s in samples]
            bystrandcheckbox = [Bonito.Checkbox(s.bystrand, style=minimal_style) for s in samples]
            byhapcheckbox = [Bonito.Checkbox(s.byhap, style=minimal_style) for s in samples]
            for (k, rtf) in enumerate(rowtextfield)
                on(rtf.value) do val
                    samples[k].row[] = parse(Int, val)
                end
            end
            for (k, rtf) in enumerate(bystrandcheckbox)
                on(rtf.value) do val

                    samples[k].bystrand = val
                    display(samples)
                    @show samples[k]
                end
            end
            for (k, rtf) in enumerate(byhapcheckbox)
                on(rtf.value) do val
                    samples[k].byhap = val
                end
            end
            sampledf = DataFrame(
                Sample=samplecheckboxtable[].Sample[inds],
                Type=[s.tracktype for s in samples],
                ByStrand=bystrandcheckbox,
                ByHap=byhapcheckbox,
                DisplayRow=rowtextfield,
            )
            rename!(sampledf, :ByStrand => "Split by strand", :ByHap => "Split by haplotype")
            sampledf
        end


        selected_display = lift(selectedsamples) do df
            if isempty(df)
                DOM.div("No samples selected")
            else
                DOM.div(Bonito.Table(df))
            end
        end


        #### figure size controls
        plotwidthinput = TextField("1700", style=minimal_style)
        trackheightinput = TextField("400", style=minimal_style)

        plotwidth = lift(plotwidthinput.value) do val
            try
                parse(Int, val)
            catch
                1700
            end
        end
        trackheight = lift(trackheightinput.value) do val
            try
                parse(Int, val)
            catch
                400
            end
        end

        mainplot = lift(chrom, start, stop) do chrom, start, stop
            if isempty(samples[]) || (chrom == "chr1" && start == 1 && stop == 1000)
                DOM.div("")
            else
                DOM.div(browserplot(chrom, start:stop, samples[], samplecheckboxtable[], genemodels=genemodels[].gene_ivs, plotwidth=plotwidth[], trackheight=trackheight[]))
            end
        end


        return DOM.div(
            style=main_style,
            DOM.h1("SMG Genome Browser"),
            DOM.div("Genome: ", genome),
            DOM.details(
                DOM.summary(DOM.span(
                    "Select samples to plot:",
                    style="font-size: 1.5em; font-weight: 600;"
                )),
                displaytable,
                open=true
            ),
            DOM.details(
                DOM.summary(DOM.span(
                    "Selected samples",
                    style="font-size: 1.5em; font-weight: 600;"
                )),
                selected_display,
            ),
            DOM.details(
                DOM.summary(DOM.span(
                    "Adaptive Sample Selection:",
                    style="font-size: 1.5em; font-weight: 600;"
                )),
                DOM.div(
                    DOM.div("Dataset: "),
                    ast_library_dropdown,
                    DOM.div("Region: "),
                    ast_region_dropdown,
                    ast_load_button,
                    style="""
                        display: flex;
                        gap: 0.5rem;
                        align-items: center;
                    """
                ),
                open=true
            ),
            DOM.details(
                DOM.summary(DOM.span(
                    "Region Navigation",
                    style="font-size: 1.5em; font-weight: 600;"
                )),
                DOM.div(
                    DOM.div("Search Location/Gene name: "),
                    locationsearch,
                    gobutton,
                    DOM.div("Current Location: "),
                    DOM.div(lift((c, s, st) -> string(c, ":", s, "-", st), chrom, start, stop)),
                    style="""
                        display: flex;
                        gap: 0.5rem;
                        align-items: center;
                    """
                ),
                open=true
            ),
            DOM.details(
                DOM.summary(DOM.span(
                    "Figure Size",
                    style="font-size: 1.5em; font-weight: 600;"
                )),
                DOM.div(
                    DOM.div("Plot width: ", plotwidthinput),
                    DOM.div("Track height: ", trackheightinput),
                    style="""
                        display: flex;
                        gap: 0.5rem;
                        align-items: center;
                    """
                ),
                open=false
            ),
            mainplot
        )
    end
end


function startserver(app, port=8080, ip="127.0.0.1")
    println("Starting server on http://$ip:$port")
    server = Bonito.Server(app, ip, port)
    server
end

# app = smgbrowser(test_load());
# close(server)
# server = startserver(app)