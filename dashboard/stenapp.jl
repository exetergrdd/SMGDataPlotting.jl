using Bonito
using WGLMakie
using SMGDataPlotting
using GenomicFeatures
using DataFrames

function loaddatasets()
    println("Loading datasets...")
    asbam_path = "/Users/ndlo201/mnt/tf/stencilling/101225_znf808_as_hac.haplotagged.fire.bam"
    wtfile_path = "/Users/ndlo201/projects/smf/fire/hap/0307_H1_60min/fire/0307_H1_60min.fire.cram"
    wgsbam_path = "/Users/ndlo201/mnt/tf/stencilling/051225_znf808/fire/051225_znf808.fire.cram"


    df = DataFrame(Study=String[], Sample=String[], MTase=String[], AdaptiveSample=Bool[], File=String[])
    push!(df, ["ZNF808", "ZNF808_KO_WGS", "EcoGII", false, wgsbam_path])
    push!(df, ["ZNF808", "ZNF808_KO_AS", "EcoGII", true, asbam_path])
    push!(df, ["H1", "H1_WGS", "EcoGII", false, wtfile_path])

    datasets = Dict("ZNF808_AS" => asbam_path, "H1WT" => wtfile_path, "ZNF808_WGS" => wgsbam_path)

    return datasets, df
end


function loadintervals()
    println("Loading intervals...")

    targets_path = "/Users/ndlo201/projects/smf/adaptive_sampling/astargets_znf808_mer11.hg19.bed"
    mer11_path = "/Users/ndlo201/resource/human/hg19/hg19.repeats.mer11.bed"

    gene_ivs = loadgenemodels("/Users/ndlo201/resource/human/hg19/gencode.v45lift37.annotation.exondata.tsv.gz")
    target_ivs = SMGDataPlotting.bedintervals(SMGDataPlotting.loadbed(targets_path))
    mer_ivs = SMGDataPlotting.bedintervals(SMGDataPlotting.loadbed(mer11_path))

    gene_dict = Dict(GenomicFeatures.metadata(iv).genename => iv for iv in gene_ivs)
    target_dict = Dict(GenomicFeatures.metadata(iv) => iv for iv in target_ivs)
    target_names = String.(sort(collect(keys(target_dict))))

    genedata = (; gene_dict, gene_ivs)
    asdata = (; target_dict, target_names)
    beddata = Dict("MER11" => mer_ivs, "Targets" => target_ivs)

    (; genedata, asdata, beddata)
end

function update_inputs!(inputs, iv)
    inputs.chrom_input.value[] = GenomicFeatures.seqname(iv)
    inputs.start_input.value[] = string(GenomicFeatures.leftposition(iv))
    inputs.stop_input.value[] = string(GenomicFeatures.rightposition(iv))
end

function gene_search_update!(query, inputs, intervals, error_msg)
    error_msg[] = "Gene '$query'."
    @show query
    if haskey(intervals.genedata.gene_dict, query)
        iv = intervals.genedata.gene_dict[query]
        update_inputs!(inputs, iv)

    else
        error_msg[] = "Gene '$query' not found."
    end
end

function runapp()


    app = App() do

        datasets, datadf = loaddatasets()
        intervals = loadintervals()
        # Initial state
        initial_chrom = "chr12"
        initial_start = 29857351
        initial_stop = 30107351


        # UI Components
        ####### define target interval #####
        genome_input = Dropdown(["hg19", "hg38"])
        chrom_input = TextField(initial_chrom)
        start_input = TextField(string(initial_start))
        stop_input = TextField(string(initial_stop))
        target_interval_inputs = (; chrom_input, start_input, stop_input)

        #### Gene and Adaptive Sampling Search ####
        gene_search = TextField("TMTC1")
        gene_go_btn = Button("Go")
        target_dropdown = Dropdown(intervals.asdata.target_names; placeholder="Select Target")

        ### Update and Zoom ###
        update_btn = Button("Update View")
        zoom_in_btn = Button("Zoom In (+)")
        zoom_out_btn = Button("Zoom Out (-)")
        reset_zoom_btn = Button("Reset Zoom") # This resets visual zoom to current inputs
        error_msg = Observable("")

        #### Event Handlers ####
        on(gene_go_btn) do click
            query = gene_search.value[]
            gene_search_update!(query, target_interval_inputs, intervals, error_msg)

        end
        on(gene_search) do text
            query = gene_search.value[]
            gene_search_update!(query, target_interval_inputs, intervals, error_msg)
        end
        on(target_dropdown.value) do target
            iv = intervals.asdata.target_dict[target]
            update_inputs!(target_interval_inputs, iv)
        end

        on(zoom_in_btn) do click
            s = parse(Int, start_input.value[])
            e = parse(Int, stop_input.value[])
            f = div(e - s, 4)
            centre = div(e + s, 2)
            new_start = s + f
            new_stop = e - f
            if new_stop - new_start < 100
                new_start = centre - 50
                new_stop = centre + 50
            end

            start_input.value[] = string(new_start)
            stop_input.value[] = string(new_stop)
        end

        on(zoom_out_btn) do click
            s = parse(Int, start_input.value[])
            e = parse(Int, stop_input.value[])
            f = div(e - s, 4)
            start_input.value[] = string(max(1, s - f))
            stop_input.value[] = string(e + f)
        end

        # on(reset_zoom_btn) do click
        #     iv = intervals.asdata.target_dict[target_dropdown.value[]]
        #     update_inputs!(target_interval_inputs, iv)
        # end

        # return DOM.div(gene_go_btn, DOM.h1(Chrom: ", chrom_input.value[]), chrom_input, start_input, stop_input, gene_search, error_msg)

        # Layout Styles
        container_style = Styles(
            "display" => "flex",
            "flex-direction" => "column",
            "gap" => "10px",
            "padding" => "20px",
            "background-color" => "#f3f4f6",
            "border-bottom" => "1px solid #e5e7eb",
            "font-family" => "sans-serif"
        )

        row_style = Styles(
            "display" => "flex",
            "gap" => "15px",
            "align-items" => "center"
        )

        input_style = Styles(
            "padding" => "5px",
            "border" => "1px solid #ccc",
            "border-radius" => "4px"
        )

        # Apply styles to inputs
        # Bonito components might accept style/class. 
        # For now we wrap them or just rely on default.
        # But for dimensions, wrapping is safe.

        # return DOM.div(
        #     container_style,

        #     # Top Row: Navigation
        #     DOM.div(
        #         row_style,
        #         DOM.span("Location: "),
        #         DOM.div(chrom_input),
        #         DOM.div(start_input),
        #         DOM.div(stop_input),
        #         update_btn
        #     ),
        #     DOM.div(
        #         row_style,
        #         DOM.span("Search Gene:"),
        #         DOM.div(gene_search),
        #         gene_go_btn, DOM.span("| Target:"),
        #         DOM.div(target_dropdown),
        #         DOM.span("| Zoom:"),
        #         zoom_in_btn,
        #         zoom_out_btn,
        #         reset_zoom_btn, DOM.span(error_msg)
        #     )
        # )

        sampletable = Bonito.Table(datadf)
        return DOM.div(
            style=container_style,

            # Top Row: Navigation
            sampletable,
            DOM.div(
                style=row_style,
                DOM.span("Genome: ", style="font-weight: bold;"),
                DOM.div(genome_input, style="width: 100px;"),
            ),
            DOM.div(
                style=row_style,
                DOM.span("Search Gene:", style="font-weight: 600;"),
                DOM.div(gene_search, style="width: 120px;"),
                gene_go_btn, DOM.span("| Target:", style="font-weight: 600; margin-left: 10px;"),
            ),
            DOM.div(
                style=row_style,
                DOM.span("Location: ", style="font-weight: bold;"),
                DOM.div(chrom_input, style="width: 100px;"),
                DOM.div(start_input, style="width: 120px;"),
                DOM.div(stop_input, style="width: 120px;"),
                update_btn,
                DOM.div(target_dropdown, style="width: 200px;"), DOM.span("| Zoom:", style="font-weight: 600; margin-left: 10px;"),
                zoom_in_btn,
                zoom_out_btn,
                reset_zoom_btn, DOM.span(error_msg, style="color: red; margin-left: 10px; font-weight: bold;")
            ),




            # Debug / Test Output
            DOM.div(
                style=row_style,
                DOM.span("Debug: Current Chrom: "),
                DOM.span(chrom_input.value)
            )
        )
    end

end

app = runapp()

port = 9386
println("Starting server on http://localhost:$port")
server = Bonito.Server(app, "127.0.0.1", port)
wait(server)
close(server)
typeof(server)
methodswith(typeof(server))