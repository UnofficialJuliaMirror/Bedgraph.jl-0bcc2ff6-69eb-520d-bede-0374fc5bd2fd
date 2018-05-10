__precompile__()

module Bedgraph
using DataFrames

export Interval

# Orthogonality.
nucleotides(n::UnitRange{Int}) = collect(n) #Note: to me it feels unreasonable to collect a range.


struct Interval
    chrom::String
    first::Int
    last::Int
    value::Real
end

function Base.:(==)(a::Interval, b::Interval)
    return a.chrom  == b.chrom &&
           a.first == b.first &&
           a.last == b.last &&
           a.value == b.value
end

mutable struct BedgraphHeader{T} #TODO: determine what and how this will be.
    data::T
end

struct BedgraphData
    header::BedgraphHeader
    intervals::Vector{Interval}
    BedgraphData(header, intervals) = new(BedgraphHeader(header), intervals) #TODO improve.
end

function Base.convert{T<:Vector{String}}(::Type{String}, header::BedgraphHeader{T}) :: String

    str = ""
    for line in header.data
        str = string(str, line, '\n')
    end

    return str
end

function Interval(data::Vector{String})
    return convert(Interval, data)
end

function Base.convert(::Type{Interval}, data::Vector{String})
    c1, c2, c3, c4 = _convertCells(data)
    return Interval(c1, c2, c3, c4)
end

function Interval(data::String)
    return convert(Interval, data)
end

function Base.convert(::Type{Interval}, str::String)
    data = _parseLine(str)
    return convert(Interval, data)
end

function Base.convert(::Type{Vector{Interval}}, chroms::Vector{String}, interval_firsts::Vector{Int}, interval_lasts::Vector{Int}, interval_values::Vector{T}) where {T<:Real}

    # Check that arrays are of equal length.
    length(chroms) == length(interval_firsts) && length(interval_lasts) == length(interval_values) && length(chroms) == length(interval_values) || error("Vectors are of unequal lengths: chroms=$(length(chroms)), interval_firsts=$(length(interval_firsts)), interval_lasts=$(length(interval_lasts)), interval_values=$(length(interval_values))")

    N = length(chroms)

    intervals = Vector{Interval}(N)

    for i in 1:N
        intervals[i] = Interval(chroms[i], interval_firsts[i], interval_lasts[i], interval_values[i])
    end

    return intervals
end

# Check if the interval data is in the four column BED format.
function isLikeInterval(line::String) :: Bool
    return  ismatch(r"^\s*\S*(?=[A-Za-z])\S*\s+(\d+)\s+(\d+)\s+(\S*\d)\s*$", line) # Note: is like a Interval.
end

function isBrowser(line::String) :: Bool
    return  ismatch(r"^browser", lowercase(line))
end

function isComment(line::String) :: Bool
    return ismatch(r"^\s*(?:#|$)", line)
end


function seekNextInterval(io) :: Void
    seekstart(io)

    pos = position(io)
    line = ""

    while !eof(io) && !isLikeInterval(line)
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

    while !eof(io) && !isLikeInterval(line) # Note: regex is used to limit the search by exiting the loop when a line matches the bedGraph interval format.
        line = readline(io)

        if contains(line, "type=bedGraph") # Note: the interval type is REQUIRED, and must be bedGraph.
            return line
        end

    end
end

function readHeader(io) :: Vector{String}
    position(io) == 0 || seekstart(io)

    header = String[]
    line = readline(io)

    while !eof(io) && !isLikeInterval(line) # TODO: seek more rebust check.
        push!(header, line)
        line = readline(io)
    end

    return header

end

function readIntervals(io) :: Vector{Interval}
    seekNextInterval(io)

    intervals = Interval[]

    while !eof(io)
        push!(intervals, Interval(readline(io)))
    end

    return intervals

end

function read(file::AbstractString, sink=DataFrame)
    # sink = Data.stream!(Source(file), sink)
    # Data.close!(sink)

    data = open(file, "r") do io
        seekNextInterval(io)
		return readdlm(io)
	end

    sink = DataFrame(chrom=data[:,1], first=data[:,2], last=data[:,3], value=data[:,4])

    return sink
end

function compress(chroms::Vector{String}, n::Vector{Int}, values::Vector{<:Real}; right_open = true, bump_back=true) :: Vector{Interval}

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

    new_intervals = Vector{Interval}()

    for (index, range) in enumerate(ranges)
        new_interval  = Interval(compressed_chroms[index], first(range), last(range), compressed_values[index])
        push!(new_intervals, new_interval)
    end

    return bump_back ? _bumpBack(new_intervals) : new_intervals

end
compress(chrom::String, n::Vector{Int}, values::Vector{T}; right_open = true, bump_back=true) where {T<:Real} = compress(fill(chrom, length(n)), n, values, right_open = right_open, bump_back = bump_back)


function compress(n::Vector{Int}, v::Vector{T}) where {T<:Real} #TODO: deprecate.

    # chrom::Vector{String} = []
    interval_firsts::Vector{Int} = []
    interval_lasts::Vector{Int} = []
    interval_values::Vector{Real} = []

    # Start inital interval.
    # push!(chrom, c[1])
    push!(interval_firsts, n[1])

    previous_value = v[1]

    state = start(v)
    while !done(v, state)
        (value, state) = next(v, state)

        # Finish current interval and start new interval if value has changed.
        if value != previous_value
            # Push interval end.
            push!(interval_lasts, n[state-1])

            # Push interval value.
            push!(interval_values, previous_value)

            # Start new interval
            # push!(chrom, c[1])

            # Push interval start.
            push!(interval_firsts, n[state-1])

            previous_value = value
        end

        if done(v, state)

            # Push final interval end.
            push!(interval_lasts, n[state-1])

            # Push final interval value.
            push!(interval_values, value)
        end
    end

    # return (chrom, first, last, value)
    return (interval_firsts, interval_lasts, interval_values)
end

compress(n, v) =  compress(nucleotides(n), v)



function expand(intervals::Vector{Interval}; right_open=true, bump_forward=true)

    #TODO: ensure intervals are sorted with no overlap.

    if bump_forward
        intervals =  _bumpForward(intervals)
    end

    total_range =_range(intervals, right_open = right_open)

    values = Vector{Float64}(length(total_range))
    chroms = Vector{String}(length(total_range))

    for interval in intervals
        values[findin(total_range, _range(interval, right_open = right_open))] = interval.value
        chroms[findin(total_range, _range(interval, right_open = right_open))] = interval.chrom
    end

    return collect(total_range), values, chroms
end

function expand(interval_firsts::Vector{Int}, interval_lasts::Vector{Int}, interval_values::Vector{T}) where {T<:Real} #TODO: deprecate.

    # Check that array are of equal length.
    if length(interval_firsts) != length(interval_lasts) || length(interval_lasts) != length(interval_values)
        error("Unequal lengths: firsts=$(length(interval_firsts)), lasts=$(length(interval_lasts)), values=$(length(interval_values))")
    end

    nucleotides = interval_firsts[1] : interval_lasts[end]
    values = zeros(length(nucleotides))

    slide = nucleotides[1] - 1

    for n = 1:length(interval_values)

        nStart = interval_firsts[n] - slide
        nEnd = interval_lasts[n] - slide

        # if left value is greater start + 1.
        if n > 1
            if interval_values[n-1] > interval_values[n]
                nStart = nStart + 1
            end
        end

        values[nStart:nEnd] = interval_values[n]
    end

    return (nucleotides, values)
end

expand(chrom::String, interval_firsts::Vector{Int}, interval_lasts::Vector{Int}, interval_values::Vector{T}; right_open=true, bump_forward=true) where {T<:Real} = expand( fill(chrom, length(interval_firsts)), interval_firsts, interval_lasts, interval_values, right_open=right_open, bump_forward=bump_forward)
expand(chroms::Vector{String}, interval_firsts::Vector{Int}, interval_lasts::Vector{Int}, interval_values::Vector{T}; right_open=true, bump_forward=true) where {T<:Real} = expand( convert(Vector{Interval}, chroms, interval_firsts, interval_lasts, interval_values), right_open=right_open, bump_forward=bump_forward)


function generateBasicHeader(chrom::String, pos_start::Int, pos_end::Int; bump_forward=true) :: Vector{String}

    if bump_forward
        pos_start += 1
        pos_end += 1
    end

    return ["browser position $chrom:$pos_start-$pos_end", "interval type=bedGraph"]
end

function generateBasicHeader(intervals::Vector{Interval}; bump_forward=true) :: Vector{String}

    chrom = intervals[1].chrom

    if bump_forward
        pos_start = intervals[1].first + 1
        pos_end = intervals[end].last + 1
    end

    return ["browser position $chrom:$pos_start-$pos_end", "interval type=bedGraph"]
end

# chrom  first  last  value
function write(chroms::Vector{String}, interval_firsts::Vector{Int}, interval_lasts::Vector{Int}, interval_values::Vector{T} ; outfile="out.bedgraph") where {T<:Real} #TODO: deprecate

    # Check that array are of equal length.
    if length(chroms) != length(interval_firsts) || length(interval_lasts) != length(interval_values) || length(chroms) != length(interval_values)
        error("Unequal lengths: chroms=$(length(chroms)), interval_firsts=$(length(interval_firsts)), interval_lasts=$(length(interval_lasts)), interval_values=$(length(interval_values))")
    end

    open(outfile, "w") do f

        for i = 1:length(chroms)

            # Write chrom.
            Base.write(f, chroms[i])
            Base.write(f, "\t")

            # Write chrom start.
            Base.write(f, string(interval_firsts[i]))
            Base.write(f, "\t")

            # Write chrom end.
            Base.write(f, string(interval_lasts[i]))
            Base.write(f, "\t")

            # Write data value.
            Base.write(f, string(interval_values[i]))
            Base.write(f, "\n")
        end

    end

end

write(c, interval_firsts, interval_lasts, interval_values; outfile="out.bedgraph" ) =  write(chrom(c), interval_firsts, interval_lasts, interval_values; outfile="out.bedgraph") #TODO: deprecate

function Base.write(io::IO, intervals::Vector{Interval}) #Note: we assume the indexes have been bumpped and the open ends are correct.
    for interval in intervals
        Base.write(io, interval, '\n')
    end
end

function Base.write(io::IO, interval::Interval)
    # delim = '\t'
    delim = ' '
    return Base.write(io, string(interval.chrom, delim, interval.first, delim, interval.last, delim, interval.value))
end

function Base.write(io::IO, header::BedgraphHeader)
    return Base.write(io, convert(String, header))
end

function Base.write(io::IO, bedgraph::BedgraphData)
    Base.write(io, bedgraph.header)
    Base.write(io, bedgraph.intervals)
end

## Internal helper functions.
function _parseLine(line::String) ::Vector{String}
    cells::Vector{String} = filter!(!isempty, split(line, r"\s"))
end

function _convertCells(cells::Vector{String})
    length(cells) == 4 || error("Poor formatting:", cells)
    return cells[1], parse(Int, cells[2]), parse(Int, cells[3]), parse(Float64, cells[4]) #TODO: parse cell 4 as a generic Real.
end


function _bump(intervals::Vector{Interval}, b::Int) :: Vector{Interval}

    new_intervals = Vector{Interval}()

    for interval in intervals
        new_interval  = Interval(interval.chrom, interval.first + b, interval.last + b, interval.value)
        push!(new_intervals, new_interval)
    end

    return new_intervals
end
_bumpForward(intervals::Vector{Interval}) = _bump(intervals, 1)
_bumpBack(intervals::Vector{Interval}) = _bump(intervals, -1)

function _range(interval::Interval; right_open=true) :: UnitRange{Int}

    pos_start = right_open ? interval.first : interval.first + 1
    pos_end = right_open ? interval.last - 1 : interval.last

    return pos_start : pos_end
end

function _range(intervals::Vector{Interval}; right_open=true) :: UnitRange{Int}

    pos_start = _range(intervals[1], right_open=right_open)[1]
    pos_end = _range(intervals[end], right_open=right_open)[end]

    return  pos_start : pos_end
end

end # module
