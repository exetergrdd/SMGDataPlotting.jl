function makedisplaytable_directrna(allmetadata, genome)


    metadata = allmetadata[allmetadata.Genome.==genome, :]

    stackcheckboxes = [Bonito.Checkbox(false) for i in 1:size(metadata, 1)]
    blockcheckboxes = [Bonito.Checkbox(false) for i in 1:size(metadata, 1)]

    dtable = deepcopy(metadata)
    dtable.StackPlot = stackcheckboxes
    dtable.BlockPlot = blockcheckboxes

    dtable
end

function smgbrowserdirectrna(metadata)


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
        genome = Bonito.Dropdown(["hg38"], style=minimal_style)
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




        ###### Sample logic

        samplecheckboxtable = lift(genome.value) do g
            empty!(samples[])
            notify(samples)
            dtable = makedisplaytable_directrna(metadata, g)
            checkboxlistener(dtable, :StackPlot, samples)
            checkboxlistener(dtable, :BlockPlot, samples)
            dtable
        end
        displaytable = lift(d -> Bonito.Table(d[!, Not([:FileClass, :File])]), samplecheckboxtable)


        selectedsamples = lift(samples) do samples

            if isempty(samples)
                return DataFrame()
            end

            inds = [s.sampleindex for s in samples]
            rowtextfield = [TextField(string(s.row[]), style=minimal_style) for s in samples]
            bystrandcheckbox = [Bonito.Checkbox(s.bystrand, style=minimal_style) for s in samples]

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
            sampledf = DataFrame(
                Sample=samplecheckboxtable[].Sample[inds],
                Type=[s.tracktype for s in samples],
                ByStrand=bystrandcheckbox,
                DisplayRow=rowtextfield,
            )
            rename!(sampledf, :ByStrand => "Split by strand")
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
                # DOM.div(browserplot(chrom, start:stop, samples[], samplecheckboxtable[], genemodels=genemodels[].gene_ivs, plotwidth=plotwidth[], trackheight=trackheight[]))
                #             smgbrowserplot(chrom, loc, displaytracks, sampletable; genemodels=nothing, datacache=Dict{DataCacheKey,Vector{T}}(), plotwidth=1800, trackheight=200,
                # datafun=SMGDataPlotting.collectmodsfire, dataplotfun=dataplot, internaltitlefun=internaltitle_stencilling)
                DOM.div(smgbrowserplot(chrom, start:stop, samples[], samplecheckboxtable[], genemodels=genemodels[].gene_ivs, plotwidth=plotwidth[], trackheight=trackheight[]),
                    datafun=SMGDataPlotting.collectdirectrnaablock, dataplotfun=dataplot_directrna, internaltitlefun=internaltitle_directRNA)
            end
        end


        return DOM.div(
            style=main_style,
            DOM.h1("SMG Direct RNA Genome Browser"),
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

# smgbrowserdirectrna(load_directrna_sampletable("test/directrna_samples.yaml"))