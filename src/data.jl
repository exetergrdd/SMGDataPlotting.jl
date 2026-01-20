struct FireRead
    lp::Int
    rp::Int
    strand::Bool
    nucs::Vector{NTuple{2, Int}}
    msps::Vector{NTuple{3, Int}}
    haplotype::UInt8
end


struct ModFireRead
    lp::Int
    rp::Int
    strand::Bool
    mods_6mA::Vector{Int}
    mods_5mC::Vector{Int}
    mods_5hmC::Vector{Int}
    nucs::Vector{NTuple{2, Int}}
    msps::Vector{NTuple{3, Int}}
    haplotype::UInt8
end



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
    sort!(reads, by = fr -> fr.lp - fr.rp);
    levels = assignlevels(reads, haplotype=byhap    )

    ### now 
    (;reads, levels)
end



function collectmodsfire(file, chrom, loc; byhap=false)
    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapModFire())

    reads = Vector{ModFireRead}()

    mods_6mA = Vector{Int}()
    mods_5mC = Vector{Int}()
    mods_5hmC = Vector{Int}()

    for record in eachintersection(reader, chrom, loc)
        validflag(record) || continue
        processread!(record, recorddata) || continue

        mods_6mA = Int[]
        mods_5mC = Int[]
        mods_5hmC = Int[]

        for mi in ModIterator(record, recorddata)
            if mi.prob > 0.9*255
                genpos = recorddata.alignmap[mi.pos]
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
    sort!(reads, by = fr -> fr.lp - fr.rp);
    levels = assignlevels(reads, haplotype=byhap    )

    ### now 
    (;reads, levels)
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

function pilereads(file::String, chrom, loc)

    @show chrom, loc

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
            pile[k - loc.start + 1] += 1
        end
    end
    pile
end

function readintervals(file, chrom, loc)
    reader = open(HTSFileReader, file)
    recorddata = StencillingData(AuxMapMod())
    reads = UnitRange{Int}[]
    # levels = IntervalLevels()
    lvls = Int[]
    for record in eachintersection(reader, chrom, loc)
        processread!(record, recorddata)
        # l = getlevel(record, recorddata, levels)
        iv = max(leftposition(record) + 1, loc.start):min(rightposition(recorddata), loc.stop)
        push!(reads, iv)

    end
    close(reader)
    reads
end