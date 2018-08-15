# Check if the track data is in the four column BED format.
function isLikeTrack(line::String) :: Bool
    return  occursin(r"^\s*\S*(?=[A-Za-z])\S*\s+(\d+)\s+(\d+)\s+(\S*\d)\s*$", line) # Note: is like a Track.
end

function isBrowser(line::String) :: Bool
    return  occursin(r"^browser", lowercase(line))
end

function isComment(line::String) :: Bool
    return occursin(r"^\s*(?:#|$)", line)
end


function seekNextTrack(io) :: Nothing
    seekstart(io)

    pos = position(io)
    line = ""

    while !eof(io) && !isLikeTrack(line)
        pos = position(io)
        line = readline(io)
    end

    seek(io, pos)

    return nothing

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

        df = DataFrame(
            chrom = Vector{String}(),
            chrom_start = Vector{Int64}(),
            chrom_end = Vector{Int64}(),
            data_value = Vector())

        while !eof(io)
            line = readline(io)
            chrom, chrom_start, chrom_end, data_value = _convertCells(_splitLine(line))

            append!(df, DataFrame(
                chrom = chrom,
                chrom_start = chrom_start,
                chrom_end = chrom_end,
                data_value = data_value))
        end

        df
	end

    # tracks = rea


    # sink = DataFrame(chrom=data[:,1], chrom_start=data[:,2], chrom_end=data[:,3], data_value=data[:,4])

    return data
end
