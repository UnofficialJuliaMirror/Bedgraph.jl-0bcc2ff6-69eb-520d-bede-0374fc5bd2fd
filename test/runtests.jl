using Bedgraph
using DataFrames
using Base.Test

@testset "Bedgraph" begin

const chroms = ["chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19"]
const interval_firsts = [49302000, 49302300, 49302600, 49302900, 49303200, 49303500, 49303800, 49304100, 49304400]
const interval_lasts = [49302300, 49302600, 49302900, 49303200, 49303500, 49303800, 49304100, 49304400, 49304700]
const interval_values = [-1.0, -0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75, 1.00]


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

const interval1 = Interval("chr19", 49302000, 49302300, -1.0)

const parameter_line_min = "track type=bedGraph"
const parameter_line = "track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20"
const parameter_line_4 = "track type=bedGraph name=track_label description=center_label"
const parameter_line_long = "track type=bedGraph name=track_label description=center_label visibility=display_mode color=r, g, b altColor=r, g, b priority=priority autoScale=on|off alwaysZero=on|off gridDefault=on|off maxHeightPixels=max:default:min graphType=bar|points viewLimits=lower:upper yLineMark=real-value yLineOnOff=on|off windowingFunction=maximum|mean|minimum smoothingWindow=off|2-16"

const file = joinpath(@__DIR__, "data.bedgraph")
const file_headerless = joinpath(@__DIR__, "data-headerless.bedgraph")

const header = [browser1, browser2, browser3, browser4, comment1, comment2, comment3, comment4, parameter_line]
const intervals = [Interval(interval1), Interval(line2), Interval(line3), Interval(line4), Interval(line5), Interval(line6), Interval(line7), Interval(line8), Interval(line9)]

@testset "I/O" begin

@test isfile(file)

# Seek test.
open(file, "r") do io
    Bedgraph.seekNextInterval(io)
    @test position(io) == 532
    @test readline(io) == line1
end

# Check things for headerless files.
open(file_headerless, "r") do io
    Bedgraph.seekNextInterval(io)
    @test position(io) == 0
    @test readline(io) == line1
end

open(file, "r") do io # Note: reading intervals first to check seek.
    @test Bedgraph.readIntervals(io) ==  intervals
    @test Bedgraph.readHeader(io) == header
end

open(file_headerless, "r") do io # Note: reading intervals first to check seek.
    @test Bedgraph.readIntervals(io) == intervals
    @test Bedgraph.readHeader(io) == []
end

# Read test.
df = Bedgraph.read(file)

@test size(df) == (9, 4)

@test df[:chrom] == chroms
@test df[:first] == interval_firsts
@test df[:last] == interval_lasts
@test df[:value] == interval_values

# Write test.
output_file = tempname() * ".bedgraph"
info(output_file)

try
    Bedgraph.write(chroms, interval_firsts, interval_lasts, interval_values, outfile=output_file)

    reloaded_df = Bedgraph.read(output_file)

    @test df == reloaded_df
finally
    rm(output_file)
end


output_file = tempname() * ".bedgraph"
info(output_file)
bedgraph = Bedgraph.BedgraphData(header, intervals)

try
    open(output_file, "w") do io
        write(io, bedgraph)
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
        write(io, Bedgraph.BedgraphData(Bedgraph.generateBasicHeader("chr19", intervals[1].first, intervals[end].last, bump_forward=false), intervals))
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

@test Bedgraph.isLikeInterval("1 2 3 4") == false
@test Bedgraph.isLikeInterval(parameter_line) == false
@test Bedgraph.isLikeInterval(parameter_line_4) == false
@test Bedgraph.isLikeInterval(parameter_line_min) == false
@test Bedgraph.isLikeInterval(parameter_line_long) == false
@test Bedgraph.isLikeInterval(line1) == true
@test Bedgraph.isLikeInterval(line1_2) == true
@test Bedgraph.isLikeInterval(line1_3) == true
@test Bedgraph.isLikeInterval(line1_4) == true
@test Bedgraph.isLikeInterval(line1_5) == true
@test Bedgraph.isLikeInterval(line1_6) == true
@test Bedgraph.isLikeInterval(line1_7) == true


@test Bedgraph.isLikeInterval(line_other_space) == true
@test Bedgraph.isLikeInterval(line_other) == true

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

@test Interval(line1) == interval1
@test convert(Interval, line1) == interval1

@test Interval(cells1) == interval1
@test convert(Interval, cells1) == interval1

@test_throws MethodError convert(Interval, String(line1, " ", "extra_cell")) #TODO: determine difference between MethodError and ErrorException.
@test_throws ErrorException convert(Interval, [cells1; "extra_cell"])

@test convert(Vector{Interval}, chroms, interval_firsts, interval_lasts, interval_values) == intervals

end #testset

@testset "Internal Helpers" begin

@test Bedgraph._range(interval1) == interval1.first : interval1.last - 1
@test Bedgraph._range(interval1, right_open=false) == (interval1.first + 1 ) : interval1.last

@test Bedgraph._range(intervals) == interval1.first : Interval(line9).last - 1
@test Bedgraph._range(intervals, right_open=false) == interval1.first + 1 : Interval(line9).last


bumped_intervals = Bedgraph._bumpForward(intervals)
@test bumped_intervals[1].first == (intervals[1].first + 1)
@test bumped_intervals[1].last == (intervals[1].last + 1)

bumped_intervals = Bedgraph._bumpBack(intervals)
@test bumped_intervals[1].first == (intervals[1].first - 1)
@test bumped_intervals[1].last == (intervals[1].last - 1)


end #testset

@testset "Utilities" begin

# Original expansion and compression test.
(n, expanded_value) = Bedgraph.expand(interval_firsts, interval_lasts, interval_values)
(compressed_first, compressed_last, compressed_value) = Bedgraph.compress(n, expanded_value)
@test interval_firsts == compressed_first
@test interval_lasts == compressed_last
@test interval_values == compressed_value

# Expansion and compression test.
n, expanded_value = Bedgraph.expand(intervals, right_open=true)
compressed_intervals = Bedgraph.compress("chr19", n, expanded_value, right_open=true)
@test compressed_intervals == intervals

# Expansion and compression of Intervals.
n, expanded_value = Bedgraph.expand(intervals, right_open=false)
compressed_intervals = Bedgraph.compress("chr19", n, expanded_value, right_open=false)
@test compressed_intervals == intervals

# Expansion and compression of Arrays via convert.
n, expanded_value = Bedgraph.expand("chr19", interval_firsts, interval_lasts, interval_values)
compressed_intervals = Bedgraph.compress("chr19", n, expanded_value)
@test compressed_intervals == intervals

end #testset

end # total testset
