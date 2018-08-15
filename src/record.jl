export Record

struct Record
    chrom::String
    chrom_start::Int
    chrom_end::Int
    value::Real
end

function Base.:(==)(a::Record, b::Record)
    return a.chrom  == b.chrom &&
           a.chrom_start == b.chrom_start &&
           a.chrom_end == b.chrom_end &&
           a.value == b.value
end

function Record(data::Vector{String})
    return convert(Record, data)
end

function Base.convert(::Type{Record}, data::Vector{String})
    c1, c2, c3, c4 = _convertCells(data)
    return Record(c1, c2, c3, c4)
end

function Record(data::String)
    return convert(Record, data)
end

function Base.convert(::Type{Record}, str::AbstractString)
    data = _splitLine(str)
    return convert(Record, data)
end

## Internal helper functions.
function _splitLine(line::String) ::Vector{String}
    cells::Vector{String} = filter!(!isempty, split(line, r"\s"))
end

function _convertCells(cells::Vector{String})
    length(cells) == 4 || error("Poor formatting:", cells)
    return cells[1], parse(Int, cells[2]), parse(Int, cells[3]), parse(Float64, cells[4]) #TODO: parse cell 4 as a generic Real.
end

function chrom(record::Record)::String
    return record.chrom
end

function chromstart(record::Record)::Int
    return record.chromstart
end

function chromend(record::Record)::Int
    return record.chromend
end

function value(record::Record)::Real
    return record.value
end
