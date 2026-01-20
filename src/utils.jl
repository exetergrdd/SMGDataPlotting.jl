### functions for fiber layout

mutable struct IntervalLevels{T, F}
    levels::Dict{T, Dict{Int, Vector{UnitRange{Int}}}}
    maxlevel::Dict{T, Int}
    levelgroup::F
end

IntervalLevels(levelgroup = (r, rd) -> 0x00) = IntervalLevels(Dict{UInt8, Dict{Int, Vector{UnitRange{Int}}}}(), Dict{UInt8, Int}(), levelgroup)

# recordlevels()    = IntervalLevels(Dict{UInt8, Dict{Int, Vector{Int}}}(), Dict{UInt8, Int}(), (r, rdata) -> 0x00)
# haplotypelevels() = IntervalLevels(Dict{UInt8, Dict{Int, Vector{Int}}}(), Dict{UInt8, Int}(), (r, rdata) -> haplotype(r, rdata, 0x00))




function assignlevels(firereads::Vector{T}; haplotype=false) where {T}
    leveldata = IntervalLevels()
    if haplotype
        [getlevel(fr.lp:fr.rp, leveldata, fr.haplotype) for fr in firereads]
    else
        [getlevel(fr.lp:fr.rp, leveldata) for fr in firereads]
    end
end


function getlevel(record, recorddata, levels::IntervalLevels)
    iv = (SMGReader.leftposition(record)+1):SMGReader.rightposition(recorddata)
    key = levels.levelgroup(record, recorddata)

   getlevel(iv, levels, key)
end

function getlevel(iv, levels, key=0x00)
    if !haskey(levels.levels, key)
        levels.levels[key] = Dict{Int, UnitRange{Int}}()
        levels.maxlevel[key] = 0
    end

    level = assign_interval_level!(iv, levels.levels[key], levels.maxlevel[key])
    levels.maxlevel[key] = max(level, levels.maxlevel[key])

    level
end
"""
    assign_interval_level!(interval, level_dict, max_level)

Assigns a genomic interval to the first available level in `level_dict` where it does not intersect with any existing intervals. 
If all levels up to `max_level` have intersections, it creates a new level.

# Arguments
- `interval`: The genomic interval to be added.
- `level_dict`: A dictionary where keys are levels and values are lists of intervals.
- `max_level`: The maximum level to check for available space.

# Returns
- The level to which the interval was assigned.

# Example
# ```julia
# level_dict = Dict{Int, Vector{UnitRange{Int}}}()
# interval = 1:1
# max_level = 1
# assigned_level = assign_interval_level!(interval, level_dict, max_level)
# ```
"""
function assign_interval_level!(interval, level_dict, max_level)
    # Iterate through each level up to max_level
    for level in 1:max_level
        if !haskey(level_dict, level)
            # If the level does not exist, create it and add the interval
            level_dict[level] = [interval]
            return level
        else
            # Check for intersections with existing intervals in the current level
            intersects = false
            for existing_interval in level_dict[level]
                if !isempty(intersect(interval, existing_interval))
                    intersects = true
                    break
                end
            end
            # If no intersections, add the interval to the current level
            if !intersects
                push!(level_dict[level], interval)
                return level
            end
        end    
    end

    # If no suitable level is found, create a new level
    level_dict[max_level + 1] = [interval]
    return max_level + 1
end

function binvector(H::Vector{T}, δ) where {T}
    n = length(H)

    w = cld(n, δ)
    AH = zeros(Float64, w)
    te = zeros(Int, w)
    for i = 1:n
        wi = cld(i, δ)
        te[wi] += 1
        AH[wi] += H[i]
    end
    if te[end] < 0.75δ
        AH[end] += AH[end-1]
        te[end] += te[end-1]
    end
    AH./te
end