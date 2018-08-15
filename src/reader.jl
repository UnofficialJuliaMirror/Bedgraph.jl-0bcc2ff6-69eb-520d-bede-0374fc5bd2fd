# Check if the record data is in the four column BED format.
function isLikeRecord(line::String) :: Bool
    return  occursin(r"^\s*\S*(?=[A-Za-z])\S*\s+(\d+)\s+(\d+)\s+(\S*\d)\s*$", line) # Note: is like a record.
end

function isBrowser(line::String) :: Bool
    return  occursin(r"^browser", lowercase(line))
end

function isComment(line::String) :: Bool
    return occursin(r"^\s*(?:#|$)", line)
end


function seekNextRecord(io) :: Nothing
    seekstart(io)

    pos = position(io)
    line = ""

    while !eof(io) && !isLikeRecord(line)
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

    while !eof(io) && !isLikeRecord(line) # Note: regex is used to limit the search by exiting the loop when a line matches the bedGraph record format.
        line = readline(io)

        if contains(line, "type=bedGraph") # Note: the track type is REQUIRED, and must be bedGraph.
            return line
        end

    end
end



function readRecords(io) :: Vector{Record}
    seekNextRecord(io)

    records = Vector{Record}()

    while !eof(io)
        push!(records, Record(readline(io)))
    end

    return records

end

function read(file::AbstractString, sink=DataFrame)
    # sink = Data.stream!(Source(file), sink)
    # Data.close!(sink)

    data = open(file, "r") do io
        seekNextRecord(io)

        df = DataFrame(
            chrom = Vector{String}(),
            first = Vector{Int64}(),
            last = Vector{Int64}(),
            value = Vector())

        while !eof(io)
            line = readline(io)
            chrom, first, last, value = _convertCells(_splitLine(line))

            append!(df, DataFrame(
                chrom = chrom,
                first = first,
                last = last,
                value = value))
        end

        df
	end

    # records = rea


    # sink = DataFrame(chrom=data[:,1], first=data[:,2], last=data[:,3], value=data[:,4])

    return data
end
