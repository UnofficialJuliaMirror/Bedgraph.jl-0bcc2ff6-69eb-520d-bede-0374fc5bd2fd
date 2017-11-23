using Bedgraph
using DataFrames
using Base.Test

@testset "Bedgraph" begin

const chrom = ["chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19"]
const chromStart = [49302000, 49302300, 49302600, 49302900, 49303200, 49303500, 49303800, 49304100, 49304400]
const chromEnd = [49302300, 49302600, 49302900, 49303200, 49303500, 49303800, 49304100, 49304400, 49304700]
const dataValue = [-1.0, -0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75, 1.00]

file = joinpath(@__DIR__, "data.bedgraph")

@test isfile(file)

open(file, "r") do io
    Bedgraph.seekNextTrack(io)
    @test position(io) == 532
    @test readline(io) == "chr19 49302000 49302300 -1.0"
end

# Read test.
df = Bedgraph.read(file)

@test size(df) == (9,4)

@test df[:chrom] == chrom
@test df[:chromStart] == chromStart
@test df[:chromEnd] == chromEnd
@test df[:dataValue] == dataValue

# Expansion and compression test.
(n, expanded_dataValue) = Bedgraph.expand(chromStart, chromEnd, dataValue)

(compressed_chromStart,compressed_chromEnd,compressed_dataValue) = Bedgraph.compress(n,expanded_dataValue)

@test chromStart == compressed_chromStart
@test chromEnd == compressed_chromEnd
@test dataValue == compressed_dataValue

# Write test.
output_file = tempname() * ".bedgraph"

try
    Bedgraph.write(chrom, chromStart, chromEnd, dataValue, outfile=output_file)

    reloaded_df = Bedgraph.read(output_file)

    @test df == reloaded_df
finally
    gc()
    rm(output_file)
end

end
