mutable struct BedgraphHeader{T} #TODO: determine what and how this will be.
    data::T
end

struct BedgraphData
    header::BedgraphHeader
    tracks::Vector{Track}
    BedgraphData(header, tracks) = new(BedgraphHeader(header), tracks) #TODO improve.
end

function Base.convert{T<:Vector{String}}(::Type{String}, header::BedgraphHeader{T}) :: String

    str = ""
    for line in header.data
        str = string(str, line, '\n')
    end

    return str
end

function generateBasicHeader(chrom::String, pos_start::Int, pos_end::Int; bump_forward=true) :: Vector{String}

    if bump_forward
        pos_start += 1
        pos_end += 1
    end

    return ["browser position $chrom:$pos_start-$pos_end", "track type=bedGraph"]
end

function generateBasicHeader(tracks::Vector{Track}; bump_forward=true) :: Vector{String}

    chrom = tracks[1].chrom

    pos_start = tracks[1].chrom_start
    pos_end = tracks[end].chrom_end

    if bump_forward
        pos_start = pos_start + 1
        pos_end = pos_end + 1
    end

    return ["browser position $chrom:$pos_start-$pos_end", "track type=bedGraph"]
end

function readHeader(io) :: Vector{String}
    position(io) == 0 || seekstart(io)

    header = String[]
    line = readline(io)

    while !eof(io) && !isLikeTrack(line) # TODO: seek more rebust check.
        push!(header, line)
        line = readline(io)
    end

    return header

end

function Base.write(io::IO, header::BedgraphHeader)
    return Base.write(io, convert(String, header))
end
