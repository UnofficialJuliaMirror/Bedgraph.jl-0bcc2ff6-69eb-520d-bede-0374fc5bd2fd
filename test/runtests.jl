using Bedgraph
using DataFrames
using Base.Test

@testset "Bedgraph" begin

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
const tracks = [Track(track1), Track(line2), Track(line3), Track(line4), Track(line5), Track(line6), Track(line7), Track(line8), Track(line9)]

@testset "I/O" begin

@test isfile(file)

# Seek test.
open(file, "r") do io
    Bedgraph.seekNextTrack(io)
    @test position(io) == 532
    @test readline(io) == line1
end

# Check things for headerless files.
open(file_headerless, "r") do io
    Bedgraph.seekNextTrack(io)
    @test position(io) == 0
    @test readline(io) == line1
end

open(file, "r") do io # Note: reading tracks first to check seek.
    @test Bedgraph.readTracks(io) ==  tracks
	@test Bedgraph._readHeader(io) == header
	@test read(io, Bedgraph.BedgraphHeader{Vector{String}}).data == header
end

open(file_headerless, "r") do io # Note: reading tracks first to check seek.
    @test Bedgraph.readTracks(io) == tracks
	@test Bedgraph._readHeader(io) == []
    @test read(io, Bedgraph.BedgraphHeader{Vector{String}}).data == []
end

# Read test.
df = Bedgraph.read(file)

@test size(df) == (9,4)

@test df[:chrom] == chroms
@test df[:chrom_start] == chrom_starts
@test df[:chrom_end] == chrom_ends
@test df[:data_value] == data_values

# Write test.
output_file = tempname() * ".bedgraph"
info(output_file)

try
    Bedgraph.write(chroms, chrom_starts, chrom_ends, data_values, outfile=output_file)

    reloaded_df = Bedgraph.read(output_file)

    @test df == reloaded_df
finally
    rm(output_file)
end


output_file = tempname() * ".bedgraph"
info(output_file)
header = Bedgraph.BedgraphHeader(Bedgraph.generateBasicHeader(tracks))

try
    open(output_file, "w") do io
        write(io, header, tracks)
    end
    # @test   readstring(file) ==  readstring(output_file) # differnces in float representation, but otherwise hold the same information.
    #TODO: explicitly test that files hold the same information.
finally
    rm(output_file)
end

output_file = tempname() * ".bedgraph"
info(output_file)

try
    open(output_file, "w") do io
        write(io, Bedgraph.BedgraphHeader(Bedgraph.generateBasicHeader("chr19", tracks[1].chrom_start, tracks[end].chrom_end, bump_forward=false)), tracks)
    end
    # @test   readstring(file) ==  readstring(output_file) # differnces in float representation, but otherwise hold the same information.
    #TODO: explicitly test that files hold the same information.
finally
    rm(output_file)
end

end #testset


@testset "Matching" begin

@test Bedgraph.isComment(comment1)
@test Bedgraph.isBrowser(browser3)

@test Bedgraph.isLikeTrack("1 2 3 4") == false
@test Bedgraph.isLikeTrack(parameter_line) == false
@test Bedgraph.isLikeTrack(parameter_line_4) == false
@test Bedgraph.isLikeTrack(parameter_line_min) == false
@test Bedgraph.isLikeTrack(parameter_line_long) == false
@test Bedgraph.isLikeTrack(line1) == true
@test Bedgraph.isLikeTrack(line1_2) == true
@test Bedgraph.isLikeTrack(line1_3) == true
@test Bedgraph.isLikeTrack(line1_4) == true
@test Bedgraph.isLikeTrack(line1_5) == true
@test Bedgraph.isLikeTrack(line1_6) == true
@test Bedgraph.isLikeTrack(line1_7) == true


@test Bedgraph.isLikeTrack(line_other_space) == true
@test Bedgraph.isLikeTrack(line_other) == true

end #testset


@testset "Parsing" begin

@test Bedgraph._parseLine(line1) == cells1
@test Bedgraph._parseLine(line1_2) == cells1
@test Bedgraph._parseLine(line1_3) == cells1
@test Bedgraph._parseLine(line1_4) == cells1
@test Bedgraph._parseLine(line1_5) == cells1
@test Bedgraph._parseLine(line1_6) == cells1
@test Bedgraph._parseLine(line1_7) == cells1

end #testset


@testset "Conversion" begin

@test_throws ErrorException Bedgraph._convertCells([cells1; "extra_cell"]) == cells1

c1, c2, c3, c4 = Bedgraph._convertCells(Bedgraph._parseLine(line1))

@test typeof(c1) == String
@test typeof(c2) <: Int
@test typeof(c3) <: Int
@test typeof(c4) <: Real

@test Track(line1) == track1
@test convert(Track, line1) == track1

@test Track(cells1) == track1
@test convert(Track, cells1) == track1

@test_throws MethodError convert(Track, String(line1, " ", "extra_cell")) #TODO: determine difference between MethodError and ErrorException.
@test_throws ErrorException convert(Track, [cells1; "extra_cell"])

@test convert(Vector{Track}, chroms, chrom_starts, chrom_ends, data_values) == tracks

end #testset

@testset "Internal Helpers" begin

@test Bedgraph._range(track1) == track1.chrom_start : track1.chrom_end - 1
@test Bedgraph._range(track1, right_open=false) == (track1.chrom_start + 1 ) : track1.chrom_end

@test Bedgraph._range(tracks) == track1.chrom_start : Track(line9).chrom_end - 1
@test Bedgraph._range(tracks, right_open=false) == track1.chrom_start + 1 : Track(line9).chrom_end


bumped_tracks = Bedgraph._bumpForward(tracks)
@test bumped_tracks[1].chrom_start == (tracks[1].chrom_start + 1)
@test bumped_tracks[1].chrom_end == (tracks[1].chrom_end + 1)

bumped_tracks = Bedgraph._bumpBack(tracks)
@test bumped_tracks[1].chrom_start == (tracks[1].chrom_start - 1)
@test bumped_tracks[1].chrom_end == (tracks[1].chrom_end - 1)


end #testset

@testset "Utilities" begin

# Original expansion and compression test.
(n, expanded_data_value) = Bedgraph.expand(chrom_starts, chrom_ends, data_values)
(compressed_chrom_start,compressed_chrom_end,compressed_data_value) = Bedgraph.compress(n,expanded_data_value)
@test chrom_starts == compressed_chrom_start
@test chrom_ends == compressed_chrom_end
@test data_values == compressed_data_value

# Expansion and compression test.
n, expanded_data_value = Bedgraph.expand(tracks, right_open=true)
compressed_tracks = Bedgraph.compress("chr19", n, expanded_data_value, right_open=true)
@test compressed_tracks == tracks

# Expansion and compression of Tracks.
n, expanded_data_value = Bedgraph.expand(tracks, right_open=false)
compressed_tracks = Bedgraph.compress("chr19", n, expanded_data_value, right_open=false)
@test compressed_tracks == tracks

# Expansion and compression of Arrays via convert.
n, expanded_data_value = Bedgraph.expand("chr19", chrom_starts, chrom_ends, data_values)
compressed_tracks = Bedgraph.compress("chr19", n, expanded_data_value)
@test compressed_tracks == tracks

end #testset

end # total testset
