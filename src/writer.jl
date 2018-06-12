# chrom  chrom_start  chrom_end  data_value
function write(chroms::Vector{String}, chrom_starts::Vector{Int}, chrom_ends::Vector{Int}, data_values::Vector{T} ; outfile="out.bedgraph") where {T<:Real} #TODO: deprecate

    # Check that array are of equal length.
    if length(chroms) != length(chrom_starts) || length(chrom_ends) != length(data_values) || length(chroms) != length(data_values)
        error("Unequal lengths: chroms=$(length(chroms)), chrom_starts=$(length(chrom_starts)), chrom_ends=$(length(chrom_ends)), data_values=$(length(data_values))")
    end

    open(outfile, "w") do f

        for i = 1:length(chroms)

            # Write chrom.
            Base.write(f, chroms[i])
            Base.write(f, "\t")

            # Write chrom start.
            Base.write(f, string(chrom_starts[i]))
            Base.write(f, "\t")

            # Write chrom end.
            Base.write(f, string(chrom_ends[i]))
            Base.write(f, "\t")

            # Write data value.
            Base.write(f, string(data_values[i]))
            Base.write(f, "\n")
        end

    end

end

write(c, chrom_starts, chrom_ends, data_values; outfile="out.bedgraph" ) =  write(chrom(c), chrom_starts, chrom_ends, data_values; outfile="out.bedgraph") #TODO: deprecate

function Base.write(io::IO, tracks::Vector{Track}) #Note: we assume the indexes have been bumpped and the open ends are correct.
    for track in tracks
        Base.write(io, track, '\n')
    end
end

function Base.write(io::IO, track::Track)
    # delim = '\t'
    delim = ' '
    return Base.write(io, string(track.chrom, delim, track.chrom_start, delim, track.chrom_end, delim, track.data_value))
end

function Base.write(io::IO, header::BedgraphHeader{Vector{String}}, tracks::Vector{Track})
    Base.write(io, header)
    Base.write(io, tracks)
end
