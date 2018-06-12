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



function Base.convert(::Type{Vector{Track}}, chroms::Vector{String}, chrom_starts::Vector{Int}, chrom_ends::Vector{Int}, data_values::Vector{T}) where {T<:Real}

    # Check that arrays are of equal length.
    length(chroms) == length(chrom_starts) && length(chrom_ends) == length(data_values) && length(chroms) == length(data_values) || error("Vectors are of unequal lengths: chroms=$(length(chroms)), chrom_starts=$(length(chrom_starts)), chrom_ends=$(length(chrom_ends)), data_values=$(length(data_values))")

    N = length(chroms)

    tracks = Vector{Track}(N)

    for i in 1:N
        tracks[i] = Track(chroms[i], chrom_starts[i], chrom_ends[i], data_values[i])
    end

    return tracks
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
    chroms = Vector{String}(length(total_range))

    for track in tracks
        values[findin(total_range, _range(track, right_open = right_open))] = track.data_value
        chroms[findin(total_range, _range(track, right_open = right_open))] = track.chrom
    end

    return collect(total_range), values, chroms
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

expand(chrom::String, chrom_starts::Vector{Int}, chrom_ends::Vector{Int}, data_values::Vector{T}; right_open=true, bump_forward=true) where {T<:Real} = expand( fill(chrom, length(chrom_starts)), chrom_starts, chrom_ends, data_values, right_open=right_open, bump_forward=bump_forward)
expand(chroms::Vector{String}, chrom_starts::Vector{Int}, chrom_ends::Vector{Int}, data_values::Vector{T}; right_open=true, bump_forward=true) where {T<:Real} = expand( convert(Vector{Track}, chroms, chrom_starts, chrom_ends, data_values), right_open=right_open, bump_forward=bump_forward)
