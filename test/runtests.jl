using Bedgraph
using DataFrames
using Test

module Bag
using Bedgraph
const chroms = ["chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19"]
const chrom_starts = [49302000, 49302300, 49302600, 49302900, 49303200, 49303500, 49303800, 49304100, 49304400]
const chrom_ends = [49302300, 49302600, 49302900, 49303200, 49303500, 49303800, 49304100, 49304400, 49304700]
const data_values = [-1.0, -0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75, 1.00]


const browser1 = "browser position chr19:49302001-49304701"
const browser2 = "browser hide all"
const browser3 = "browser pack refGene encodeRegions"
const browser4 = "browser full altGraph"

const comment1 = "#	300 base wide bar graph, autoScale is on by default == graphing"
const comment2 = "#	limits will dynamically change to always show full range of data"
const comment3 = "#	in viewing window, priority = 20 positions this as the second graph"
const comment4 = "#	Note, zero-relative, half-open coordinate system in use for bedGraph format"

# space separated lines.
const line1 = "chr19 49302000 49302300 -1.0"
const line2 = "chr19 49302300 49302600 -0.75"
const line3 = "chr19 49302600 49302900 -0.50"
const line4 = "chr19 49302900 49303200 -0.25"
const line5 = "chr19 49303200 49303500 0.0"
const line6 = "chr19 49303500 49303800 0.25"
const line7 = "chr19 49303800 49304100 0.50"
const line8 = "chr19 49304100 49304400 0.75"
const line9 = "chr19 49304400 49304700 1.00"

const line_other_space = "2R 8225773 8226043 -0.426032509896305"
const line_other = "2R	8225773	8226043	-0.426032509896305"

# Varaiations of line 1.
const line1_2 = "chr19   49302000    49302300    -1.0" # tab separated.
const line1_3 = "chr19  49302000     49302300        -1.0" # mix of tabs and spaces.
const line1_4 = " chr19 49302000 49302300 -1.0" # space at start.
const line1_5 = "    chr19 49302000 49302300 -1.0" # tab at start.
const line1_6 = "chr19 49302000 49302300 -1.0 " # space at end.
const line1_7 = "chr19 49302000 49302300 -1.0    " # tab at end.

const cells1 = ["chr19", "49302000", "49302300", "-1.0"]

const track1 = Track("chr19", 49302000, 49302300, -1.0)

const parameter_line_min = "track type=bedGraph"
const parameter_line = "track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20"
const parameter_line_4 = "track type=bedGraph name=track_label description=center_label"
const parameter_line_long = "track type=bedGraph name=track_label description=center_label visibility=display_mode color=r,g,b altColor=r,g,b priority=priority autoScale=on|off alwaysZero=on|off gridDefault=on|off maxHeightPixels=max:default:min graphType=bar|points viewLimits=lower:upper yLineMark=real-value yLineOnOff=on|off windowingFunction=maximum|mean|minimum smoothingWindow=off|2-16"

const file = joinpath(@__DIR__, "data.bedgraph")
const file_headerless = joinpath(@__DIR__, "data-headerless.bedgraph")

const header = [browser1, browser2, browser3, browser4, comment1, comment2, comment3, comment4, parameter_line]
const tracks = [track1, Track(line2), Track(line3), Track(line4), Track(line5), Track(line6), Track(line7), Track(line8), Track(line9)]

end # module bag

@testset "Bedgraph" begin

@testset "I/O" begin

@test isfile(Bag.file)
@test isfile(Bag.file_headerless)

# Seek test.
open(Bag.file, "r") do io
    Bedgraph.seekNextTrack(io)
    @test position(io) == 532
    @test readline(io) == Bag.line1
end

# Check things for headerless Bag.files.
open(Bag.file_headerless, "r") do io
    Bedgraph.seekNextTrack(io)
    @test position(io) == 0
    @test readline(io) == Bag.line1
end

open(Bag.file, "r") do io # Note: reading tracks first to check seek.
    @test Bedgraph.readTracks(io) ==  Bag.tracks
	@test Bedgraph._readHeader(io) == Bag.header
	@test read(io, Bedgraph.BedgraphHeader{Vector{String}}).data == Bag.header
end

open(Bag.file_headerless, "r") do io # Note: reading tracks first to check seek.
    @test Bedgraph.readTracks(io) == Bag.tracks
	@test Bedgraph._readHeader(io) == []
    @test read(io, Bedgraph.BedgraphHeader{Vector{String}}).data == []
end

# Read test.
df = Bedgraph.read(Bag.file)

@test size(df) == (9,4)

@test df[:chrom] == Bag.chroms
@test df[:chrom_start] == Bag.chrom_starts
@test df[:chrom_end] == Bag.chrom_ends
@test df[:data_value] == Bag.data_values

# Write test.
outputfile1 = tempname() * ".bedgraph"
@info "outputfile:" outputfile1

try
    Bedgraph.write(Bag.chroms, Bag.chrom_starts, Bag.chrom_ends, Bag.data_values, outfile=outputfile1)

    reloaded_df = Bedgraph.read(outputfile1)

    @test df == reloaded_df
finally
    rm(outputfile1)
end


outputfile2 = tempname() * ".bedgraph"
@info "outputfile:" outputfile2
header = Bedgraph.BedgraphHeader(Bedgraph.generateBasicHeader(Bag.tracks))

try
    open(outputfile2, "w") do io
        write(io, header, Bag.tracks)
    end
    # @test   readstring(Bag.file) ==  readstring(outputfile) # differnces in float representation, but otherwise hold the same information.
    #TODO: explicitly test that Bag.files hold the same information.
finally
    rm(outputfile2)
end

outputfile3 = tempname() * ".bedgraph"
@info "outputfile:" outputfile3

try
    open(outputfile3, "w") do io
        header = Bedgraph.BedgraphHeader(Bedgraph.generateBasicHeader("chr19", Bag.tracks[1].chrom_start, Bag.tracks[end].chrom_end, bump_forward=false))
        write(io, header, Bag.tracks)
    end
    # @test   readstring(Bag.file) ==  readstring(outputfile) # differnces in float representation, but otherwise hold the same information.
    #TODO: explicitly test that Bag.files hold the same information.
finally
    rm(outputfile3)
end

end #testset I/O


@testset "Matching" begin

@test Bedgraph.isComment(Bag.comment1)
@test Bedgraph.isBrowser(Bag.browser3)

@test Bedgraph.isLikeTrack("1 2 3 4") == false
@test Bedgraph.isLikeTrack(Bag.parameter_line) == false
@test Bedgraph.isLikeTrack(Bag.parameter_line_4) == false
@test Bedgraph.isLikeTrack(Bag.parameter_line_min) == false
@test Bedgraph.isLikeTrack(Bag.parameter_line_long) == false
@test Bedgraph.isLikeTrack(Bag.line1) == true
@test Bedgraph.isLikeTrack(Bag.line1_2) == true
@test Bedgraph.isLikeTrack(Bag.line1_3) == true
@test Bedgraph.isLikeTrack(Bag.line1_4) == true
@test Bedgraph.isLikeTrack(Bag.line1_5) == true
@test Bedgraph.isLikeTrack(Bag.line1_6) == true
@test Bedgraph.isLikeTrack(Bag.line1_7) == true


@test Bedgraph.isLikeTrack(Bag.line_other_space) == true
@test Bedgraph.isLikeTrack(Bag.line_other) == true

end #testset Matching


@testset "Parsing" begin

@test Bedgraph._splitLine(Bag.line1) == Bag.cells1
@test Bedgraph._splitLine(Bag.line1_2) == Bag.cells1
@test Bedgraph._splitLine(Bag.line1_3) == Bag.cells1
@test Bedgraph._splitLine(Bag.line1_4) == Bag.cells1
@test Bedgraph._splitLine(Bag.line1_5) == Bag.cells1
@test Bedgraph._splitLine(Bag.line1_6) == Bag.cells1
@test Bedgraph._splitLine(Bag.line1_7) == Bag.cells1

end #testset Parsing


@testset "Conversion" begin

@test_throws ErrorException Bedgraph._convertCells([Bag.cells1; "extra_cell"]) == Bag.cells1

c1, c2, c3, c4 = Bedgraph._convertCells(Bedgraph._splitLine(Bag.line1))

@test typeof(c1) == String
@test typeof(c2) <: Int
@test typeof(c3) <: Int
@test typeof(c4) <: Real

@test Track(Bag.line1) == Bag.track1
@test convert(Track, Bag.line1) == Bag.track1

@test Track(Bag.cells1) == Bag.track1
@test convert(Track, Bag.cells1) == Bag.track1

@test_throws MethodError convert(Track, String(Bag.line1, " ", "extra_cell")) #TODO: determine difference between MethodError and ErrorException.
@test_throws ErrorException convert(Track, [Bag.cells1; "extra_cell"])

@test convert(Vector{Track}, Bag.chroms, Bag.chrom_starts, Bag.chrom_ends, Bag.data_values) == Bag.tracks

end #testset Conversion

@testset "Internal Helpers" begin

@test Bedgraph._range(Bag.track1) == Bag.track1.chrom_start : Bag.track1.chrom_end - 1
@test Bedgraph._range(Bag.track1, right_open=false) == (Bag.track1.chrom_start + 1 ) : Bag.track1.chrom_end

@test Bedgraph._range(Bag.tracks) == Bag.track1.chrom_start : Track(Bag.line9).chrom_end - 1
@test Bedgraph._range(Bag.tracks, right_open=false) == Bag.track1.chrom_start + 1 : Track(Bag.line9).chrom_end


bumped_tracks = Bedgraph._bumpForward(Bag.tracks)
@test bumped_tracks[1].chrom_start == (Bag.tracks[1].chrom_start + 1)
@test bumped_tracks[1].chrom_end == (Bag.tracks[1].chrom_end + 1)

bumped_tracks = Bedgraph._bumpBack(Bag.tracks)
@test bumped_tracks[1].chrom_start == (Bag.tracks[1].chrom_start - 1)
@test bumped_tracks[1].chrom_end == (Bag.tracks[1].chrom_end - 1)


end #testset Internal Helpers

@testset "Utilities" begin

# Original expansion and compression test.
(n, expanded_data_value) = Bedgraph.expand(Bag.chrom_starts, Bag.chrom_ends, Bag.data_values)
(compressed_chrom_start,compressed_chrom_end,compressed_data_value) = Bedgraph.compress(n,expanded_data_value)
@test Bag.chrom_starts == compressed_chrom_start
@test Bag.chrom_ends == compressed_chrom_end
@test Bag.data_values == compressed_data_value

# Expansion and compression test.
n, expanded_data_value = Bedgraph.expand(Bag.tracks, right_open=true)
compressed_tracks = Bedgraph.compress("chr19", n, expanded_data_value, right_open=true)
@test compressed_tracks == Bag.tracks

# Expansion and compression of Tracks.
n, expanded_data_value = Bedgraph.expand(Bag.tracks, right_open=false)
compressed_tracks = Bedgraph.compress("chr19", n, expanded_data_value, right_open=false)
@test compressed_tracks == Bag.tracks

# Expansion and compression of Arrays via convert.
n, expanded_data_value = Bedgraph.expand("chr19", Bag.chrom_starts, Bag.chrom_ends, Bag.data_values)
compressed_tracks = Bedgraph.compress("chr19", n, expanded_data_value)
@test compressed_tracks == Bag.tracks

end #testset Utilities

end # total testset
