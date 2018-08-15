function _bump(records::Vector{Record}, b::Int) :: Vector{Record}

    new_records = Vector{Record}()

    for record in records
        new_record  = Record(record.chrom, record.chrom_start + b, record.chrom_end + b, record.data_value)
        push!(new_records, new_record)
    end

    return new_records
end
_bumpForward(records::Vector{Record}) = _bump(records, 1)
_bumpBack(records::Vector{Record}) = _bump(records, -1)

function _range(record::Record; right_open=true) :: UnitRange{Int}

    pos_start = right_open ? record.chrom_start : record.chrom_start + 1
    pos_end = right_open ? record.chrom_end - 1 : record.chrom_end

    return pos_start : pos_end
end

function _range(records::Vector{Record}; right_open=true) :: UnitRange{Int}

    pos_start = _range(records[1], right_open=right_open)[1]
    pos_end = _range(records[end], right_open=right_open)[end]

    return  pos_start : pos_end
end



function Base.convert(::Type{Vector{Record}}, chroms::Vector{String}, chrom_starts::Vector{Int}, chrom_ends::Vector{Int}, data_values::Vector{T}) where {T<:Real}

    # Check that arrays are of equal length.
    length(chroms) == length(chrom_starts) && length(chrom_ends) == length(data_values) && length(chroms) == length(data_values) || error("Vectors are of unequal lengths: chroms=$(length(chroms)), chrom_starts=$(length(chrom_starts)), chrom_ends=$(length(chrom_ends)), data_values=$(length(data_values))")

    records = Vector{Record}()

    for (chrom, chrom_start, chrom_end, data_value) in zip(chroms, chrom_starts, chrom_ends, data_values)
        push!(records,  Record(chrom, chrom_start, chrom_end, data_value))
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
    chrom_starts::Vector{Int} = []
    chrom_ends::Vector{Int} = []
    data_values::Vector{Real} = []

    # Start inital record.
    # push!(chrom, c[1])
    push!(chrom_starts, n[1])

    previous_value = v[1]

    next = iterate(v)
    while next !== nothing
        (value, state) = next

        # Finish current record and start new record if value has changed.
        if value != previous_value
            # Push record end.
            push!(chrom_ends, n[state-1])

            # Push record value.
            push!(data_values, previous_value)

            # Start new record
            # push!(chrom, c[1])

            # Push record start.
            push!(chrom_starts, n[state-1])

            previous_value = value
        end

        next = iterate(v, state)

        if next == nothing
            # Push final record end.
            push!(chrom_ends, n[state-1])

            # Push final record value.
            push!(data_values, value)
        end
    end

    # return (chrom, chrom_start, chrom_end, data_value)
    return (chrom_starts, chrom_ends, data_values)
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
        values[indexin(_range(record, right_open = right_open), total_range)] .= record.data_value
        chroms[indexin(_range(record, right_open = right_open), total_range)] .= record.chrom
    end

    return collect(total_range), values, chroms
end

function expand(chrom_starts::Vector{Int}, chrom_ends::Vector{Int}, data_values::Vector{T}) where {T<:Real} #TODO: deprecate.

    # Check that array are of equal length.
    if length(chrom_starts) != length(chrom_ends) || length(chrom_ends) != length(data_values)
        error("Unequal lengths: chrom_starts=$(length(chrom_starts)), chrom_ends=$(length(chrom_ends)), data_values=$(length(data_values))")
    end

    nucleotides = chrom_starts[1] : chrom_ends[end]
    values = zeros(T, length(nucleotides))

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

        values[nStart : nEnd] .= data_values[n]
    end

    return (nucleotides, values)
end

expand(chrom::String, chrom_starts::Vector{Int}, chrom_ends::Vector{Int}, data_values::Vector{T}; right_open=true, bump_forward=true) where {T<:Real} = expand( fill(chrom, length(chrom_starts)), chrom_starts, chrom_ends, data_values, right_open=right_open, bump_forward=bump_forward)
expand(chroms::Vector{String}, chrom_starts::Vector{Int}, chrom_ends::Vector{Int}, data_values::Vector{T}; right_open=true, bump_forward=true) where {T<:Real} = expand( convert(Vector{Record}, chroms, chrom_starts, chrom_ends, data_values), right_open=right_open, bump_forward=bump_forward)
