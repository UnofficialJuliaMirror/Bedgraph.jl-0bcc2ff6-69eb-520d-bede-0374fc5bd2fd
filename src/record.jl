export Track

mutable struct Track
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

## Internal helper functions.
function _parseLine(line::String) ::Vector{String}
    cells::Vector{String} = filter!(!isempty, split(line, r"\s"))
end

function _convertCells(cells::Vector{String})
    length(cells) == 4 || error("Poor formatting:", cells)
    return cells[1], parse(Int, cells[2]), parse(Int, cells[3]), parse(Float64, cells[4]) #TODO: parse cell 4 as a generic Real.
end
