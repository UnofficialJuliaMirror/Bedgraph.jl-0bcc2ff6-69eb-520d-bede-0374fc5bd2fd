module Bedgraph
using DataFrames
# using DataArrays

# Orthogonality.
chrom(c::Array{String,1}) = c
chrom(c::SubString{String}) = string(c)
chrom(c::Array{Any,1}) = convert(Array{String,1}, c)
chrom(c::DataArrays.DataArray{Any,1}) = chrom(dropna(c))
nucleotides(n::Array{Int,1}) = n
nucleotides(n::UnitRange{Int}) = collect(n) #Note: to me it feels unreasonable to collect a range.
nucleotides(n::DataArrays.DataArray{Any,1}) = convert(Array{Int,1}, n)
dataValues(v::DataArrays.DataArray{Any,1}) = convert(Array{Float64,1}, v)

function read(file::AbstractString, sink=DataFrame)
    # sink = Data.stream!(Source(file), sink)
    # Data.close!(sink)

    data = readdlm(file)

    sink = DataFrame(chrom=data[:,1], chromStart=data[:,2], chromEnd=data[:,3], dataValue=data[:,4])

    return sink
end

function compress(n::Array{Int,1}, v::Array{Float64,1})

    # chrom::Array{String,1} = []
    chromStart::Array{Int,1} = []
    chromEnd::Array{Int,1} = []
    dataValue::Array{Float64,1} = []

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

function expand(chromStart::Array{Int,1}, chromEnd::Array{Int,1}, dataValue::Array{Float64,1})

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
function write(chrom::Array{String,1}, chromStart::Array{Int,1}, chromEnd::Array{Int,1}, dataValue::Array{Float64,1} ; outfile="out.bedgraph")

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

end # module
