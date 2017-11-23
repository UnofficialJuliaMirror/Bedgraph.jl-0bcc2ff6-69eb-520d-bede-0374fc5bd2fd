using Bedgraph
using DataFrames
using Base.Test

@testset "Bedgraph" begin

const chrom = ["chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19"]
const chromStart = [49302000, 49302300, 49302600, 49302900, 49303200, 49303500, 49303800, 49304100, 49304400]
const chromEnd = [49302300, 49302600, 49302900, 49303200, 49303500, 49303800, 49304100, 49304400, 49304700]
const dataValue = [-1.0, -0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75, 1.00]

const comment = "#	300 base wide bar graph, autoScale is on by default == graphing"
const browser = "browser position chr19:49302001-49304701"

const line1 = "chr19 49302000 49302300 -1.0" # space separated.
const line2 = "chr19   49302000    49302300    -1.0" # tab separated.
const line3 = "chr19  49302000     49302300        -1.0" # mix of tabs and spaces.
const line4 = " chr19 49302000 49302300 -1.0" # space at start.
const line5 = "    chr19 49302000 49302300 -1.0" # tab at start.
const line6 = "chr19 49302000 49302300 -1.0 " # space at end.
const line7 = "chr19 49302000 49302300 -1.0    " # tab at end.

const cells1 = ["chr19", "49302000", "49302300", "-1.0"] :: Vector{String}

const parameter_line = "track type=bedGraph"
const parameter_line_short = "track type=bedGraph name=track_label description=center_label"
const parameter_line_long = "track type=bedGraph name=track_label description=center_label visibility=display_mode color=r,g,b altColor=r,g,b priority=priority autoScale=on|off alwaysZero=on|off gridDefault=on|off maxHeightPixels=max:default:min graphType=bar|points viewLimits=lower:upper yLineMark=real-value yLineOnOff=on|off windowingFunction=maximum|mean|minimum smoothingWindow=off|2-16"

const file = joinpath(@__DIR__, "data.bedgraph")


@testset "I/O" begin

@test isfile(file)

# Seek test.
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

end #testset


@testset "Matching" begin

@test Bedgraph.isComment(comment)
@test Bedgraph.isBrowser(browser)

@test Bedgraph.isLikeTrack("1 2 3 4") == false
@test Bedgraph.isLikeTrack(parameter_line) == false
@test Bedgraph.isLikeTrack(parameter_line_short) == false
@test Bedgraph.isLikeTrack(parameter_line_long) == false
@test Bedgraph.isLikeTrack(line1) == true
@test Bedgraph.isLikeTrack(line2) == true
@test Bedgraph.isLikeTrack(line3) == true
@test Bedgraph.isLikeTrack(line4) == true
@test Bedgraph.isLikeTrack(line5) == true
@test Bedgraph.isLikeTrack(line6) == true
@test Bedgraph.isLikeTrack(line7) == true

end #testset


@testset "Parsing" begin

@test Bedgraph.parseLine(line1) == cells1
@test Bedgraph.parseLine(line2) == cells1
@test Bedgraph.parseLine(line3) == cells1
@test Bedgraph.parseLine(line4) == cells1
@test Bedgraph.parseLine(line5) == cells1
@test Bedgraph.parseLine(line6) == cells1
@test Bedgraph.parseLine(line7) == cells1


end #testset


@testset "Utilities" begin

# Expansion and compression test.
(n, expanded_dataValue) = Bedgraph.expand(chromStart, chromEnd, dataValue)

(compressed_chromStart,compressed_chromEnd,compressed_dataValue) = Bedgraph.compress(n,expanded_dataValue)

@test chromStart == compressed_chromStart
@test chromEnd == compressed_chromEnd
@test dataValue == compressed_dataValue

end #testset


end # total testset
