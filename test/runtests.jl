using Bedgraph
using DataFrames
using Test

module Bag
using Bedgraph
const chroms = ["chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19"]
const firsts = [49302000, 49302300, 49302600, 49302900, 49303200, 49303500, 49303800, 49304100, 49304400]
const chrom_ends = [49302300, 49302600, 49302900, 49303200, 49303500, 49303800, 49304100, 49304400, 49304700]
const values = [-1.0, -0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75, 1.00]


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

const record1 = Record("chr19", 49302000, 49302300, -1.0)

const parameter_line_min = "track type=bedGraph"
const parameter_line = "track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20"
const parameter_line_4 = "track type=bedGraph name=track_label description=center_label"
const parameter_line_long = "track type=bedGraph name=track_label description=center_label visibility=display_mode color=r,g,b altColor=r,g,b priority=priority autoScale=on|off alwaysZero=on|off gridDefault=on|off maxHeightPixels=max:default:min graphType=bar|points viewLimits=lower:upper yLineMark=real-value yLineOnOff=on|off windowingFunction=maximum|mean|minimum smoothingWindow=off|2-16"

const file = joinpath(@__DIR__, "data.bedgraph")
const file_headerless = joinpath(@__DIR__, "data-headerless.bedgraph")

const header = [browser1, browser2, browser3, browser4, comment1, comment2, comment3, comment4, parameter_line]
const records = [record1, Record(line2), Record(line3), Record(line4), Record(line5), Record(line6), Record(line7), Record(line8), Record(line9)]

end # module bag

@testset "Bedgraph" begin

@testset "I/O" begin

@test isfile(Bag.file)
@test isfile(Bag.file_headerless)

# Seek test.
open(Bag.file, "r") do io
    Bedgraph.seekNextRecord(io)
    @test position(io) == 532
    @test readline(io) == Bag.line1
end

# Check things for headerless Bag.files.
open(Bag.file_headerless, "r") do io
    Bedgraph.seekNextRecord(io)
    @test position(io) == 0
    @test readline(io) == Bag.line1
end

open(Bag.file, "r") do io # Note: reading records first to check seek.
    @test Bedgraph.readRecords(io) ==  Bag.records
	@test Bedgraph._readHeader(io) == Bag.header
	@test read(io, Bedgraph.BedgraphHeader{Vector{String}}).data == Bag.header
end

open(Bag.file_headerless, "r") do io # Note: reading records first to check seek.
    @test Bedgraph.readRecords(io) == Bag.records
	@test Bedgraph._readHeader(io) == []
    @test read(io, Bedgraph.BedgraphHeader{Vector{String}}).data == []
end

# Read test.
df = Bedgraph.read(Bag.file)

@test size(df) == (9,4)

@test df[:chrom] == Bag.chroms
@test df[:first] == Bag.firsts
@test df[:chrom_end] == Bag.chrom_ends
@test df[:value] == Bag.values

# Write test.
outputfile1 = tempname() * ".bedgraph"
@info "outputfile:" outputfile1

try
    Bedgraph.write(Bag.chroms, Bag.firsts, Bag.chrom_ends, Bag.values, outfile=outputfile1)

    reloaded_df = Bedgraph.read(outputfile1)

    @test df == reloaded_df
finally
    rm(outputfile1)
end


outputfile2 = tempname() * ".bedgraph"
@info "outputfile:" outputfile2
header = Bedgraph.BedgraphHeader(Bedgraph.generateBasicHeader(Bag.records))

try
    open(outputfile2, "w") do io
        write(io, header, Bag.records)
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
        header = Bedgraph.BedgraphHeader(Bedgraph.generateBasicHeader("chr19", Bag.records[1].first, Bag.records[end].chrom_end, bump_forward=false))
        write(io, header, Bag.records)
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

@test Bedgraph.isLikeRecord("1 2 3 4") == false
@test Bedgraph.isLikeRecord(Bag.parameter_line) == false
@test Bedgraph.isLikeRecord(Bag.parameter_line_4) == false
@test Bedgraph.isLikeRecord(Bag.parameter_line_min) == false
@test Bedgraph.isLikeRecord(Bag.parameter_line_long) == false
@test Bedgraph.isLikeRecord(Bag.line1) == true
@test Bedgraph.isLikeRecord(Bag.line1_2) == true
@test Bedgraph.isLikeRecord(Bag.line1_3) == true
@test Bedgraph.isLikeRecord(Bag.line1_4) == true
@test Bedgraph.isLikeRecord(Bag.line1_5) == true
@test Bedgraph.isLikeRecord(Bag.line1_6) == true
@test Bedgraph.isLikeRecord(Bag.line1_7) == true


@test Bedgraph.isLikeRecord(Bag.line_other_space) == true
@test Bedgraph.isLikeRecord(Bag.line_other) == true

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

@test Record(Bag.line1) == Bag.record1
@test convert(Record, Bag.line1) == Bag.record1

@test Record(Bag.cells1) == Bag.record1
@test convert(Record, Bag.cells1) == Bag.record1

@test_throws MethodError convert(Record, String(Bag.line1, " ", "extra_cell")) #TODO: determine difference between MethodError and ErrorException.
@test_throws ErrorException convert(Record, [Bag.cells1; "extra_cell"])

@test convert(Vector{Record}, Bag.chroms, Bag.firsts, Bag.chrom_ends, Bag.values) == Bag.records

end #testset Conversion

@testset "Internal Helpers" begin

@test Bedgraph._range(Bag.record1) == Bag.record1.first : Bag.record1.chrom_end - 1
@test Bedgraph._range(Bag.record1, right_open=false) == (Bag.record1.first + 1 ) : Bag.record1.chrom_end

@test Bedgraph._range(Bag.records) == Bag.record1.first : Record(Bag.line9).chrom_end - 1
@test Bedgraph._range(Bag.records, right_open=false) == Bag.record1.first + 1 : Record(Bag.line9).chrom_end


bumped_records = Bedgraph._bumpForward(Bag.records)
@test bumped_records[1].first == (Bag.records[1].first + 1)
@test bumped_records[1].chrom_end == (Bag.records[1].chrom_end + 1)

bumped_records = Bedgraph._bumpBack(Bag.records)
@test bumped_records[1].first == (Bag.records[1].first - 1)
@test bumped_records[1].chrom_end == (Bag.records[1].chrom_end - 1)


end #testset Internal Helpers

@testset "Utilities" begin

# Original expansion and compression test.
(n, expanded_value) = Bedgraph.expand(Bag.firsts, Bag.chrom_ends, Bag.values)
(compressed_first,compressed_chrom_end,compressed_value) = Bedgraph.compress(n,expanded_value)
@test Bag.firsts == compressed_first
@test Bag.chrom_ends == compressed_chrom_end
@test Bag.values == compressed_value

# Expansion and compression test.
n, expanded_value = Bedgraph.expand(Bag.records, right_open=true)
compressed_records = Bedgraph.compress("chr19", n, expanded_value, right_open=true)
@test compressed_records == Bag.records

# Expansion and compression of Records.
n, expanded_value = Bedgraph.expand(Bag.records, right_open=false)
compressed_records = Bedgraph.compress("chr19", n, expanded_value, right_open=false)
@test compressed_records == Bag.records

# Expansion and compression of Arrays via convert.
n, expanded_value = Bedgraph.expand("chr19", Bag.firsts, Bag.chrom_ends, Bag.values)
compressed_records = Bedgraph.compress("chr19", n, expanded_value)
@test compressed_records == Bag.records

end #testset Utilities

end # total testset
