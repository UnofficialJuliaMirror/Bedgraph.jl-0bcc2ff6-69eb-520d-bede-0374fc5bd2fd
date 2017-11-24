module Bedgraph
using DataFrames
# using DataArrays

export Track

# Orthogonality.
chrom(c::Vector{String}) = c
chrom(c::SubString{String}) = string(c)
chrom(c::Vector{Any}) = convert(Vector{String}, c)
chrom(c::DataArrays.DataArray{Any,1}) = chrom(dropna(c))
nucleotides(n::Vector{Int}) = n
nucleotides(n::UnitRange{Int}) = collect(n) #Note: to me it feels unreasonable to collect a range.
nucleotides(n::DataArrays.DataArray{Any,1}) = convert(Vector{Int}, n)
dataValues(v::DataArrays.DataArray{Any,1}) = convert(Vector{T}, v) where {T<:Real}


struct Track
    chrom::String
    chrom_start::Int
    chrom_end::Int
    data_value::Real
end

function Base.:(==)(a::Track, b::Track)
    return a.chrom  == b.chrom &&
           a.chrom_start == b.chrom_start &&
           a.chrom_end == b.chrom_end &&
           a.data_value == b.data_value
end

function Track(data::Vector{String})
    return convert(Track, data)
end

function Base.convert(::Type{Track}, data::Vector{String})
    c1, c2, c3, c4 = _convertCells(data)
    return Track(c1, c2, c3, c4)
end

function Track(data::String)
    return convert(Track, data)
end

function Base.convert(::Type{Track}, str::String)
    data = _parseLine(str)
    return convert(Track, data)
end

# Check if the track data in four column BED format.
function isLikeTrack(line::String) :: Bool
    return  ismatch(r"^\s*([A-Za-z]+\S*)\s+(\d+)\s+(\d+)\s+(\S*\d)\s*$", line) # Note: is like a Track.
end

function isBrowser(line::String) :: Bool
    return  ismatch(r"^browser", lowercase(line))
end

function isComment(line::String) :: Bool
    return ismatch(r"^\s*(?:#|$)", line)
end

#
function seekNextTrack(io) :: Void
    seekstart(io)

    pos = position(io)
    line = ""

    while !eof(io) && !isLikeTrack(line)
        pos = position(io)
        line = readline(io)
    end

    seek(io, pos)

    nothing

end

# Note: all options are placed in a single line separated by spaces.
function readParameters(io) :: String
    seekstart(io)

    pos = position(io)

    while !eof(io) && !isLikeTrack(line) # Note: regex is used to limit the search by exiting the loop when a line matches the bedGraph track format.
        line = readline(io)

        if contains(line, "type=bedGraph") # Note: the track type is REQUIRED, and must be bedGraph.
            return line
        end

    end
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

function readTracks(io) :: Vector{Track}
    seekNextTrack(io)

    tracks = Track[]

    while !eof(io)
        push!(tracks, Track(readline(io)))
    end

    return tracks

end

function read(file::AbstractString, sink=DataFrame)
    # sink = Data.stream!(Source(file), sink)
    # Data.close!(sink)

    data = open(file, "r") do io
        seekNextTrack(io)
		return readdlm(io)
	end

    sink = DataFrame(chrom=data[:,1], chromStart=data[:,2], chromEnd=data[:,3], dataValue=data[:,4])

    return sink
end



function compress(n::Vector{Int}, v::Vector{T}) where {T<:Real}

    # chrom::Vector{String} = []
    chromStart::Vector{Int} = []
    chromEnd::Vector{Int} = []
    dataValue::Vector{Real} = []

    # Start inital track.
    # push!(chrom, c[1])
    push!(chromStart, n[1])

    previous_value = v[1]

    state = start(v)
    while !done(v, state)
        (value, state) = next(v, state)

        # Finish current track and start new track if value has changed.
        if value != previous_value
            # Push track end.
            push!(chromEnd, n[state-1])

            # Push track value.
            push!(dataValue, previous_value)

            # Start new track
            # push!(chrom, c[1])

            # Push track start.
            push!(chromStart, n[state-1])

            previous_value = value
        end

        if done(v, state)

            # Push final track end.
            push!(chromEnd, n[state-1])

            # Push final track value.
            push!(dataValue, value)
        end
    end

    # return (chrom, chromStart, chromEnd, dataValue)
    return (chromStart, chromEnd, dataValue)
end

compress(n, v) =  compress(nucleotides(n), v)

function expand(chromStart::Vector{Int}, chromEnd::Vector{Int}, dataValue::Vector{T}) where {T<:Real}

    # Check that array are of equal length.
    if length(chromStart) != length(chromEnd) || length(chromEnd) != length(dataValue)
        error("Unequal lengths: chromStart=$(length(chromStart)), chromEnd=$(length(chromEnd)), dataValue=$(length(dataValue))")
    end

    nucleotides = chromStart[1] : chromEnd[end]
    values = zeros(length(nucleotides))

    slide = nucleotides[1] - 1

    for n = 1:length(dataValue)

        nStart = chromStart[n] - slide
        nEnd = chromEnd[n] - slide

        # if left value is greater start + 1.
        if n > 1
            if dataValue[n-1] > dataValue[n]
                nStart = nStart + 1
            end
        end

        values[nStart:nEnd] = dataValue[n]
    end

    return (nucleotides, values)
end

expand(chromStart::DataArrays.DataArray{Any,1}, chromEnd::DataArrays.DataArray{Any,1}, dataValue::DataArrays.DataArray{Any,1}) = expand(nucleotides(chromStart), nucleotides(chromEnd), dataValues(dataValue))

# chrom  chromStart  chromEnd  dataValue
function write(chrom::Vector{String}, chromStart::Vector{Int}, chromEnd::Vector{Int}, dataValue::Vector{T} ; outfile="out.bedgraph") where {T<:Real}

    # Check that array are of equal length.
    if length(chrom) != length(chromStart) || length(chromEnd) != length(dataValue) || length(chrom) != length(dataValue)
        error("Unequal lengths: chrom=$(length(chrom)), chromStart=$(length(chromStart)), chromEnd=$(length(chromEnd)), dataValue=$(length(dataValue))")
    end

    open(outfile, "w") do f

        for i = 1:length(chrom)

            # Write chrom.
            Base.write(f, chrom[i])
            Base.write(f, "\t")

            # Write chrom start.
            Base.write(f, string(chromStart[i]))
            Base.write(f, "\t")

            # Write chrom end.
            Base.write(f, string(chromEnd[i]))
            Base.write(f, "\t")

            # Write data value.
            Base.write(f, string(dataValue[i]))
            Base.write(f, "\n")
        end

    end

end

write(c, chromStart, chromEnd, dataValue; outfile="out.bedgraph" ) =  write(chrom(c), chromStart, chromEnd, dataValue; outfile="out.bedgraph")

## Internal helper functions.
function _parseLine(line::String) ::Vector{String}
    cells::Vector{String} = filter!(!isempty, split(line, r"\s"))
end

function _convertCells(cells::Vector{String})
    length(cells) == 4 || error("Poor line formatting:", cells)
    return cells[1], parse(Int, cells[2]), parse(Int, cells[3]), parse(Float64, cells[4]) #TODO: parse cell 4 as a generic Real.
end

function _range(track::Track; right_open=true) :: UnitRange{Int}

    pos_start = right_open ? track.chrom_start : track.chrom_start + 1
    pos_end = right_open ? track.chrom_end - 1 : track.chrom_end

    return pos_start : pos_end
end

function _range(tracks::Vector{Track}; right_open=true) :: UnitRange{Int}

    pos_start = _range(tracks[1], right_open=right_open)[1]
    pos_end = _range(tracks[end], right_open=right_open)[end]

    return  pos_start : pos_end
end

end # module
