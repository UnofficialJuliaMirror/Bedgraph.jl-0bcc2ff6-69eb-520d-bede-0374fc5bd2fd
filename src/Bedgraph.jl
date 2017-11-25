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
data_values(v::DataArrays.DataArray{Any,1}) = convert(Vector{T}, v) where {T<:Real}


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

function Base.convert(::Type{Vector{Track}}, chroms::Vector{String}, chrom_starts::Vector{Int}, chrom_ends::Vector{Int}, data_values::Vector{T}) where {T<:Real}

    # Check that arrays are of equal length.
    length(chroms) == length(chrom_starts) && length(chrom_ends) == length(data_values) && length(chroms) == length(data_values) || error("Unequal lengths: chroms=$(length(chroms)), chrom_starts=$(length(chrom_starts)), chrom_ends=$(length(chrom_ends)), data_values=$(length(data_values))")

    N = length(chroms)

    tracks = Vector{Track}(N)

    for i in 1:N
        tracks[i] = Track(chroms[i], chrom_starts[i], chrom_ends[i], data_values[i])
    end

    return tracks
end

# Check if the track data is in the four column BED format.
function isLikeTrack(line::String) :: Bool
    return  ismatch(r"^\s*([A-Za-z]+\S*)\s+(\d+)\s+(\d+)\s+(\S*\d)\s*$", line) # Note: is like a Track.
end

function isBrowser(line::String) :: Bool
    return  ismatch(r"^browser", lowercase(line))
end

function isComment(line::String) :: Bool
    return ismatch(r"^\s*(?:#|$)", line)
end


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

    sink = DataFrame(chrom=data[:,1], chrom_start=data[:,2], chrom_end=data[:,3], data_value=data[:,4])

    return sink
end

function compress(chroms::Vector{String}, n::Vector{Int}, values::Vector{<:Real}; right_open = true, bump_back=true) :: Vector{Track}

    ranges = Vector{UnitRange{Int}}()
    compressed_values = Vector{Float64}()
    compressed_chroms = Vector{String}()

    range_start = 1
    push!(compressed_values, values[1])

    for (index, value ) in enumerate(values)
        if value != compressed_values[end]
            push!(ranges, n[range_start] : n[index - 1] )
            push!(compressed_values, value)
            push!(compressed_chroms, chroms[index])
            range_start = index
        end

        if index == length(values)
            push!(ranges, n[range_start] : n[index] )
            push!(compressed_values, value)
            push!(compressed_chroms, chroms[index])
        end
    end

    if right_open
        for (index, value) in enumerate(ranges)
            ranges[index] = first(value) : last(value) + 1
        end
    else
        for (index, value) in enumerate(ranges)
            ranges[index] = first(value) -1 : last(value)
        end
    end

    new_tracks = Vector{Track}()

    for (index, range) in enumerate(ranges)
        new_track  = Track(compressed_chroms[index], first(range), last(range), compressed_values[index])
        push!(new_tracks, new_track)
    end

    return bump_back ? _bumpBack(new_tracks) : new_tracks

end
compress(chrom::String, n::Vector{Int}, values::Vector{T}; right_open = true, bump_back=true) where {T<:Real} = compress(fill(chrom, length(n)), n, values, right_open = right_open, bump_back = bump_back)


function compress(n::Vector{Int}, v::Vector{T}) where {T<:Real} #TODO: deprecate.

    # chrom::Vector{String} = []
    chrom_starts::Vector{Int} = []
    chrom_ends::Vector{Int} = []
    data_values::Vector{Real} = []

    # Start inital track.
    # push!(chrom, c[1])
    push!(chrom_starts, n[1])

    previous_value = v[1]

    state = start(v)
    while !done(v, state)
        (value, state) = next(v, state)

        # Finish current track and start new track if value has changed.
        if value != previous_value
            # Push track end.
            push!(chrom_ends, n[state-1])

            # Push track value.
            push!(data_values, previous_value)

            # Start new track
            # push!(chrom, c[1])

            # Push track start.
            push!(chrom_starts, n[state-1])

            previous_value = value
        end

        if done(v, state)

            # Push final track end.
            push!(chrom_ends, n[state-1])

            # Push final track value.
            push!(data_values, value)
        end
    end

    # return (chrom, chrom_start, chrom_end, data_value)
    return (chrom_starts, chrom_ends, data_values)
end

compress(n, v) =  compress(nucleotides(n), v)



function expand(tracks::Vector{Track}; right_open=true, bump_forward=true)

    #TODO: ensure tracks are sorted with no overlap.

    if bump_forward
        tracks =  _bumpForward(tracks)
    end

    total_range =_range(tracks, right_open = right_open)

    values = Vector{Float64}(length(total_range))

    for track in tracks
        values[findin(total_range, _range(track, right_open = right_open))] = track.data_value
    end

    return collect(total_range), values
end

function expand(chrom_starts::Vector{Int}, chrom_ends::Vector{Int}, data_values::Vector{T}) where {T<:Real} #TODO: deprecate.

    # Check that array are of equal length.
    if length(chrom_starts) != length(chrom_ends) || length(chrom_ends) != length(data_values)
        error("Unequal lengths: chrom_starts=$(length(chrom_starts)), chrom_ends=$(length(chrom_ends)), data_values=$(length(data_values))")
    end

    nucleotides = chrom_starts[1] : chrom_ends[end]
    values = zeros(length(nucleotides))

    slide = nucleotides[1] - 1

    for n = 1:length(data_values)

        nStart = chrom_starts[n] - slide
        nEnd = chrom_ends[n] - slide

        # if left value is greater start + 1.
        if n > 1
            if data_values[n-1] > data_values[n]
                nStart = nStart + 1
            end
        end

        values[nStart:nEnd] = data_values[n]
    end

    return (nucleotides, values)
end

expand(chrom_starts::DataArrays.DataArray{Any,1}, chrom_ends::DataArrays.DataArray{Any,1}, data_values::DataArrays.DataArray{Any,1}) = expand(nucleotides(chrom_starts), nucleotides(chrom_ends), data_values(data_values))
expand(chrom::String, chrom_starts::Vector{Int}, chrom_ends::Vector{Int}, data_values::Vector{T}; right_open=true, bump_forward=true) where {T<:Real} = expand( fill(chrom, length(chrom_starts)), chrom_starts, chrom_ends, data_values, right_open=right_open, bump_forward=bump_forward)
expand(chroms::Vector{String}, chrom_starts::Vector{Int}, chrom_ends::Vector{Int}, data_values::Vector{T}; right_open=true, bump_forward=true) where {T<:Real} = expand( convert(Vector{Track}, chroms, chrom_starts, chrom_ends, data_values), right_open=right_open, bump_forward=bump_forward)


function generateBasicHeader(chrom::String, pos_start::Int, pos_end::Int; bump_forward=true) :: Vector{String}

    if bump_forward
        pos_start += 1
        pos_end += 1
    end

    return ["browser position $chrom:$pos_start-$pos_end", "track type=bedGraph"]
end

# chrom  chrom_start  chrom_end  data_value
function write(chroms::Vector{String}, chrom_starts::Vector{Int}, chrom_ends::Vector{Int}, data_values::Vector{T} ; outfile="out.bedgraph") where {T<:Real}

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

write(c, chrom_starts, chrom_ends, data_values; outfile="out.bedgraph" ) =  write(chrom(c), chrom_starts, chrom_ends, data_values; outfile="out.bedgraph")

## Internal helper functions.
function _parseLine(line::String) ::Vector{String}
    cells::Vector{String} = filter!(!isempty, split(line, r"\s"))
end

function _convertCells(cells::Vector{String})
    length(cells) == 4 || error("Poor line formatting:", cells)
    return cells[1], parse(Int, cells[2]), parse(Int, cells[3]), parse(Float64, cells[4]) #TODO: parse cell 4 as a generic Real.
end


function _bump(tracks::Vector{Track}, b::Int) :: Vector{Track}

    new_tracks = Vector{Track}()

    for track in tracks
        new_track  = Track(track.chrom, track.chrom_start + b, track.chrom_end + b, track.data_value)
        push!(new_tracks, new_track)
    end

    return new_tracks
end
_bumpForward(tracks::Vector{Track}) = _bump(tracks, 1)
_bumpBack(tracks::Vector{Track}) = _bump(tracks, -1)

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
