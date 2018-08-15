__precompile__()

module Bedgraph
using DataFrames
using TextParse
# using GenomicFeatures

include("record.jl")
include("helpers.jl")
include("header.jl")
include("reader.jl")
include("writer.jl")

end # module
