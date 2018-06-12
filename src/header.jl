mutable struct BedgraphHeader{T} #TODO: determine what and how this will be.
    data::T
end

function Base.convert{T<:Vector{String}}(::Type{String}, header::BedgraphHeader{T}) :: String

    str = ""
    for line in header.data
        str = string(str, line, '\n')
    end

    return str
end

function Base.convert(::Type{BedgraphHeader{Vector{String}}}, header::Vector{String})
    return BedgraphHeader(header)
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

generateBasicHeader(chrom::String, pos_start::Int, pos_end::Int; bump_forward=true) = generateBasicHeader([Track(chrom, pos_start, pos_end,0)], bump_forward=bump_forward)

function _readHeader(io) :: Vector{String}
    position(io) == 0 || seekstart(io)

    header = String[]
    line = readline(io)

    while !eof(io) && !isLikeTrack(line) # TODO: seek more rebust check.
        push!(header, line)
        line = readline(io)
    end

    return header

end

function Base.read{T}(io::IO, ::Type{BedgraphHeader{T}})
    return BedgraphHeader(_readHeader(io))
end


function Base.write(io::IO, header::BedgraphHeader{Vector{String}})
    return Base.write(io, convert(String, header))
end
