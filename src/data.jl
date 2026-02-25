"""
    FireRead

Struct representing a single read with FIRE (Fiber-seq Inferred Regulatory Elements) information.

# Fields
- `lp::Int`: Left position of the read.
- `rp::Int`: Right position of the read.
- `strand::Bool`: Strand of the read (`true` for positive, `false` for negative).
- `nucs::Vector{NTuple{2, Int}}`: Vector of nucleosome positions (start, stop).
- `msps::Vector{NTuple{3, Int}}`: Vector of MSPs (start, stop, quality).
- `haplotype::UInt8`: Haplotype assignment.
"""
struct FireRead
    lp::Int
    rp::Int
    strand::Bool
    nucs::Vector{NTuple{2,Int}}
    msps::Vector{NTuple{3,Int}}
    haplotype::UInt8
end


"""
    ModFireRead

Struct representing a single read with modified bases and FIRE information.

# Fields
- `lp::Int`: Left position of the read.
- `rp::Int`: Right position of the read.
- `strand::Bool`: Strand of the read.
- `mods_6mA::Vector{Int}`: Positions of 6mA modifications.
- `mods_5mC::Vector{Int}`: Positions of 5mC modifications.
- `mods_5hmC::Vector{Int}`: Positions of 5hmC modifications.
- `nucs::Vector{NTuple{2, Int}}`: Nucleosome positions.
- `msps::Vector{NTuple{3, Int}}`: MSP positions.
- `haplotype::UInt8`: Haplotype assignment.
"""
struct ModFireRead
    lp::Int
    rp::Int
    strand::Bool
    mods_6mA::Vector{Int}
    mods_5mC::Vector{Int}
    mods_5hmC::Vector{Int}
    nucs::Vector{NTuple{2,Int}}
    msps::Vector{NTuple{3,Int}}
    haplotype::UInt8
end

struct ReadLevels{T,L}
    reads::Vector{T}
    levels::Vector{L}
end

"""
    collectfire(file, chrom, loc; byhap=false)

Collects FIRE elements (nucleosomes and MSPs) from a BAM/CRAM file for a specific genomic region.

# Arguments
- `file`: Path to the BAM/CRAM file.
- `chrom`: Chromosome name.
- `loc`: UnitRange representing the genomic locus (e.g., `1000:2000`).
- `byhap`: Boolean to group reads by haplotype (default `false`).

# Returns
A NamedTuple `(;reads, levels)` containing the collected `FireRead` objects and their assigned levels for plotting.
"""
function collectfire(file, chrom, loc; byhap=false)
    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapModFire())

    reads = Vector{FireRead}()

    for record in eachintersection(reader, chrom, loc)
        validflag(record) || continue
        processread!(record, recorddata) || continue
        firegenomecoords = firegenome(record, recorddata)

        nucs = firegenomecoords(firenucs(record, recorddata))
        msps = firegenomecoords(firemsps(record, recorddata))

        fr = FireRead(SMGReader.leftposition(record), SMGReader.rightposition(recorddata), ispositive(record), collect(nucs), collect(msps), haplotype(record, recorddata, 0x00))
        push!(reads, fr)

    end
    close(reader)
    @time sort!(reads, by=fr -> fr.lp - fr.rp)
    # levels = assignlevels(reads, haplotype=byhap)

    # ReadLevels(reads, levels)
    reads
end



"""
    collectmodsfire(file, chrom, loc; byhap=false)

Collects modified bases (6mA, 5mC, 5hmC) and FIRE elements from a BAM/CRAM file for a specific genomic region.

# Arguments
- `file`: Path to the BAM/CRAM file.
- `chrom`: Chromosome name.
- `loc`: UnitRange representing the genomic locus.
- `byhap`: Boolean to group reads by haplotype (default `false`).

# Returns
A NamedTuple `(;reads, levels)` containing the collected `ModFireRead` objects and their assigned levels.
"""
function collectmodsfire(file, chrom, loc; byhap=false, bystrand=false, strictintersect=0)
    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapModFire())

    reads = Vector{ModFireRead}()

    mods_6mA = Vector{Int}()
    mods_5mC = Vector{Int}()
    mods_5hmC = Vector{Int}()


    for record in eachintersection(reader, chrom, loc)
        validflag(record) || continue
        processread!(record, recorddata) || continue

        if strictintersect > 0
            lp = SMGReader.leftposition(record) + 1
            rp = SMGReader.rightposition(recorddata)

            ((lp + strictintersect < first(loc)) && (last(loc) < rp + strictintersect)) || continue
        end



        mods_6mA = Int[]
        mods_5mC = Int[]
        mods_5hmC = Int[]

        for mi in ModIterator(record, recorddata)
            if mi.prob > 0.9 * 255

                genpos = genomecoords(mi.pos, record, recorddata, onebased=true)
                if !iszero(genpos)
                    if mi.mod == mod_6mA
                        push!(mods_6mA, genpos)
                    elseif mi.mod == mod_5mC
                        push!(mods_5mC, genpos)
                    elseif mi.mod == mod_5hmC
                        push!(mods_5hmC, genpos)
                    end
                end
            end
        end

        firegenomecoords = firegenome(record, recorddata)

        nucs = firegenomecoords(firenucs(record, recorddata))
        msps = firegenomecoords(firemsps(record, recorddata))

        fr = ModFireRead(SMGReader.leftposition(record), SMGReader.rightposition(recorddata), ispositive(record), mods_6mA, mods_5mC, mods_5hmC, collect(nucs), collect(msps), haplotype(record, recorddata, 0x00))
        push!(reads, fr)

    end
    close(reader)
    @time sort!(reads, by=fr -> fr.lp - fr.rp)
    # levels = assignlevels(reads, haplotype=byhap, strand=bystrand)

    ### now 
    # (; reads, levels)
    reads
end

function filterreads(readlevels, filtfun=x -> true, byhap=false)
    reads = filter(filtun, readlevels.reads)
    sort!(reads, by=fr -> fr.lp - fr.rp)
    levels = assignlevels(reads, haplotype=byhap)
    (; reads, levels)
end

function combinereads(readlevelsA, readlevelsB; byhap=false)
    reads = [readlevelsA.reads; readlevelsB.reads]
    sort!(reads, by=fr -> fr.lp - fr.rp)
    levels = assignlevels(reads, haplotype=byhap)
    (; reads, levels)
end



# function fireelements(file, chrom, loc)
#     reader = open(HTSFileReader, file)
#     recorddata = StencillingData(AuxMapModFire())

#     levels = IntervalLevels()

#     nucdata = Vector{NTuple{3, Int}}()
#     mspdata = Vector{NTuple{4, Int}}()

#     for record in eachintersection(reader, chrom, loc)
#         validflag(record) || continue
#         processread!(record, recorddata) || continue
#         firegenomecoords = firegenome(record, recorddata)
#         l = getlevel(record, recorddata, levels)

#         nucs = firegenomecoords(firenucs(record, recorddata))
#         for (start, stop) in nucs
#             push!(nucdata, (l, start, stop))
#         end

#         msps = firegenomecoords(firemsps(record, recorddata))
#         for (start, stop, qual) in msps
#             push!(mspdata, (l, start, stop, qual))
#         end
#     end

#     close(reader)

#     (; nucdata, mspdata)
# end

"""
    pilereads(file::String, chrom, loc)
    pilereads(reader, chrom, loc, recorddata)

Calculates the coverage (pileup) of reads over a specific genomic region.

# Arguments
- `file`: Path to the BAM/CRAM file.
- `chrom`: Chromosome name.
- `loc`: UnitRange representing the genomic locus.

# Returns
A `Vector{Float64}` representing the read count at each position in `loc`.
"""
function pilereads(file::String, chrom, loc)

    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapMod())
    pile = pilereads(reader, chrom, loc, recorddata)

    close(reader)
    pile

end

function pilereads(reader, chrom, loc, recorddata)
    n = length(loc)
    pile = zeros(n)
    for record in eachintersection(reader, chrom, loc)
        processread!(record, recorddata)

        iv = max(SMGReader.leftposition(record) + 1, loc.start):min(SMGReader.rightposition(recorddata), loc.stop)

        for k in iv
            pile[k-loc.start+1] += 1
        end
    end
    pile
end

"""
    readintervals(file, chrom, loc)

Extracts the genomic intervals covered by reads in a specific region.

# Arguments
- `file`: Path to the BAM/CRAM file.
- `chrom`: Chromosome name.
- `loc`: UnitRange representing the genomic locus.

# Returns
A `Vector{UnitRange{Int}}` of the intervals covered by each read.
"""
function readintervals(file, chrom, loc)
    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapMod())
    reads = UnitRange{Int}[]
    # levels = IntervalLevels()
    # lvls = Int[]
    for record in eachintersection(reader, chrom, loc)
        processread!(record, recorddata)
        # l = getlevel(record, recorddata, levels)
        iv = max(leftposition(record) + 1, loc.start):min(rightposition(recorddata), loc.stop)
        push!(reads, iv)

    end
    close(reader)
    reads
end