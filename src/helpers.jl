function _bump(records::Vector{Record}, b::Int) :: Vector{Record}

    new_records = Vector{Record}()

    for record in records
        new_record  = Record(record.chrom, record.first + b, record.last + b, record.value)
        push!(new_records, new_record)
    end

    return new_records
end
_bumpForward(records::Vector{Record}) = _bump(records, 1)
_bumpBack(records::Vector{Record}) = _bump(records, -1)

function _range(record::Record; right_open=true) :: UnitRange{Int}

    pos_start = right_open ? record.first : record.first + 1
    pos_end = right_open ? record.last - 1 : record.last

    return pos_start : pos_end
end

function _range(records::Vector{Record}; right_open=true) :: UnitRange{Int}

    pos_start = _range(records[1], right_open=right_open)[1]
    pos_end = _range(records[end], right_open=right_open)[end]

    return  pos_start : pos_end
end



function Base.convert(::Type{Vector{Record}}, chroms::Vector{String}, firsts::Vector{Int}, lasts::Vector{Int}, values::Vector{T}) where {T<:Real}

    # Check that arrays are of equal length.
    length(chroms) == length(firsts) && length(lasts) == length(values) && length(chroms) == length(values) || error("Vectors are of unequal lengths: chroms=$(length(chroms)), firsts=$(length(firsts)), lasts=$(length(lasts)), values=$(length(values))")

    records = Vector{Record}()

    for (chrom, first, last, value) in zip(chroms, firsts, lasts, values)
        push!(records,  Record(chrom, first, last, value))
    end

    return records
end


function compress(chroms::Vector{String}, n::Vector{Int}, values::Vector{<:Real}; right_open = true, bump_back=true) :: Vector{Record}

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

    new_records = Vector{Record}()

    for (index, range) in enumerate(ranges)
        new_record  = Record(compressed_chroms[index], first(range), last(range), compressed_values[index])
        push!(new_records, new_record)
    end

    return bump_back ? _bumpBack(new_records) : new_records

end
compress(chrom::String, n::Vector{Int}, values::Vector{T}; right_open = true, bump_back=true) where {T<:Real} = compress(fill(chrom, length(n)), n, values, right_open = right_open, bump_back = bump_back)


function compress(n::Vector{Int}, v::Vector{T}) where {T<:Real} #TODO: deprecate.

    # chrom::Vector{String} = []
    firsts::Vector{Int} = []
    lasts::Vector{Int} = []
    values::Vector{Real} = []

    # Start inital record.
    # push!(chrom, c[1])
    push!(firsts, n[1])

    previous_value = v[1]

    next = iterate(v)
    while next !== nothing
        (value, state) = next

        # Finish current record and start new record if value has changed.
        if value != previous_value
            # Push record end.
            push!(lasts, n[state-1])

            # Push record value.
            push!(values, previous_value)

            # Start new record
            # push!(chrom, c[1])

            # Push record start.
            push!(firsts, n[state-1])

            previous_value = value
        end

        next = iterate(v, state)

        if next == nothing
            # Push final record end.
            push!(lasts, n[state-1])

            # Push final record value.
            push!(values, value)
        end
    end

    # return (chrom, first, last, value)
    return (firsts, lasts, values)
end

compress(n, v) =  compress(collect(n), v)



function expand(records::Vector{Record}; right_open=true, bump_forward=true)

    #TODO: ensure records are sorted with no overlap.

    if bump_forward
        records =  _bumpForward(records)
    end

    total_range =_range(records, right_open = right_open)

    values = Vector{Float64}(undef, length(total_range))
    chroms = Vector{String}(undef, length(total_range))

    for record in records
        values[indexin(_range(record, right_open = right_open), total_range)] .= record.value
        chroms[indexin(_range(record, right_open = right_open), total_range)] .= record.chrom
    end

    return collect(total_range), values, chroms
end

function expand(firsts::Vector{Int}, lasts::Vector{Int}, values::Vector{T}) where {T<:Real} #TODO: deprecate.

    # Check that array are of equal length.
    if length(firsts) != length(lasts) || length(lasts) != length(values)
        error("Unequal lengths: firsts=$(length(firsts)), lasts=$(length(lasts)), values=$(length(values))")
    end

    nucleotides = firsts[1] : lasts[end]
    new_values = zeros(T, length(nucleotides))

    slide = nucleotides[1] - 1

    for n = 1:length(values)

        nStart = firsts[n] - slide
        nEnd = lasts[n] - slide

        # if left value is greater start + 1.
        if n > 1
            if values[n-1] > values[n]
                nStart = nStart + 1
            end
        end

        new_values[nStart : nEnd] .= values[n]
    end

    return (nucleotides, new_values)
end

expand(chrom::String, firsts::Vector{Int}, lasts::Vector{Int}, values::Vector{T}; right_open=true, bump_forward=true) where {T<:Real} = expand( fill(chrom, length(firsts)), firsts, lasts, values, right_open=right_open, bump_forward=bump_forward)
expand(chroms::Vector{String}, firsts::Vector{Int}, lasts::Vector{Int}, values::Vector{T}; right_open=true, bump_forward=true) where {T<:Real} = expand( convert(Vector{Record}, chroms, firsts, lasts, values), right_open=right_open, bump_forward=bump_forward)
