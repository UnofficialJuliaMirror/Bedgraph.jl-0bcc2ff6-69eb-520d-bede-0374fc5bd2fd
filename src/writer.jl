# chrom  first  last  value
function write(chroms::Vector{String}, firsts::Vector{Int}, lasts::Vector{Int}, values::Vector{T} ; outfile="out.bedgraph") where {T<:Real} #TODO: deprecate

    # Check that array are of equal length.
    if length(chroms) != length(firsts) || length(lasts) != length(values) || length(chroms) != length(values)
        error("Unequal lengths: chroms=$(length(chroms)), firsts=$(length(firsts)), lasts=$(length(lasts)), values=$(length(values))")
    end

    open(outfile, "w") do f

        for i = 1:length(chroms)

            # Write chrom.
            Base.write(f, chroms[i])
            Base.write(f, "\t")

            # Write chrom start.
            Base.write(f, string(firsts[i]))
            Base.write(f, "\t")

            # Write chrom end.
            Base.write(f, string(lasts[i]))
            Base.write(f, "\t")

            # Write data value.
            Base.write(f, string(values[i]))
            Base.write(f, "\n")
        end

    end

end

write(c, firsts, lasts, values; outfile="out.bedgraph" ) =  write(chrom(c), firsts, lasts, values; outfile="out.bedgraph") #TODO: deprecate

function Base.write(io::IO, records::Vector{Record}) #Note: we assume the indexes have been bumpped and the open ends are correct.
    for record in records
        Base.write(io, record, '\n')
    end
end

function Base.write(io::IO, record::Record)
    # delim = '\t'
    delim = ' '
    return Base.write(io, string(record.chrom, delim, record.first, delim, record.last, delim, record.value))
end

function Base.write(io::IO, header::BedgraphHeader{Vector{String}}, records::Vector{Record})
    Base.write(io, header)
    Base.write(io, records)
end
