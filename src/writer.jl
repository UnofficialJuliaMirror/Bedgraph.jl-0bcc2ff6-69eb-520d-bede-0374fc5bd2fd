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
