# Bedgraph.jl

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/CiaranOMara/Bedgraph.jl.svg?branch=master)](https://travis-ci.org/CiaranOMara/Bedgraph.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/jny2ep4u3cmly8pj/branch/master?svg=true)](https://ci.appveyor.com/project/CiaranOMara/Bedgraph-jl/branch/master)
[![Bedgraph](http://pkg.julialang.org/badges/Bedgraph_0.6.svg)](http://pkg.julialang.org/?pkg=Bedgraph)
[![codecov.io](http://codecov.io/github/CiaranOMara/Bedgraph.jl/coverage.svg?branch=master)](http://codecov.io/github/CiaranOMara/Bedgraph.jl?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/CiaranOMara/Bedgraph.jl/badge.svg?branch=master)](https://coveralls.io/github/CiaranOMara/Bedgraph.jl?branch=master)

> This project will try to follow the [semver](http://semver.org) pro forma.

## Overview
This package provides read and write support for [Bedgraph files](https://genome.ucsc.edu/goldenPath/help/bedgraph.html), as well as other useful utilities.

> **Note:**  this package does not currently handle bedGraph meta data such as the track definition or browser lines.

## Installation
Use Pkg.add("Bedgraph") in Julia to install Bedgraph.jl and its dependencies.

## Usage

### Reading and writing bedGraph files
> See source for optional `bump_back`, `bump_forward`, and `right_open` key values. These options are included in the pertinent read/write functions to handle quirks of the zero-based and half-open nature of the bedGraph format.

#### Read a bedGraph file into a DataFrame
Bedgraph.jl currently returns read data as a DataFrame.

```julia
using Bedgraph, DataFrames

df = Bedgraph.read("data.bedgraph")
```

#### Read header/meta
```julia
using Bedgraph

header = Vector{String}()
open(file, "r") do io
    header = Bedgraph.readHeader(io)
end
```

#### Read records

```julia
using Bedgraph

records = Vector{Record}()
open(file, "r") do io
    records = Bedgraph.readRecords(io)
end
```

#### Write a bedGraph file
Bedgraph.jl currently provides two options. Either vectors or a BedgraphData type can be supplied to its write function.

```julia
using Bedgraph

const chroms = ["chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19", "chr19"]
const chrom_starts = [49302000, 49302300, 49302600, 49302900, 49303200, 49303500, 49303800, 49304100, 49304400]
const chrom_ends = [49302300, 49302600, 49302900, 49303200, 49303500, 49303800, 49304100, 49304400, 49304700]
const data_values = [-1.0, -0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75, 1.00]

Bedgraph.write(chroms, chrom_starts, chrom_ends, data_values, outfile="data.bedgraph")
```


```julia
using Bedgraph

records = [Record("chr19", 49302000, 49302300, -1.0), Record("chr19", 49302300, 49302600, -1.75)]

open(output_file, "w") do io
    write(io, Bedgraph.BedgraphData(Bedgraph.generateBasicHeader("chr19", records[1].chrom_start, records[end].chrom_end, bump_forward=false), records))
end

```
### Expansion and compression of data

#### Compress data values
Compress data to chromosome coordinates of the zero-based, half-open format.

```julia
using Bedgraph

n = 49302000:49304700
expanded_data_values = [-1.0,-1.0,-1.0, ..., 1.00, 1.00, 1.00]

(compressed_chrom_starts,compressed_chrom_ends,compressed_data_values) = Bedgraph.compress(n,expanded_data_values)
```

```julia
using Bedgraph

const records = [Record("chr19", 49302000, 49302300, -1.0), Record("chr19", 49302300, 49302600, -1.75)]

compressed_records = Bedgraph.compress("chr19", n, expanded_data_value)
```

#### Expand record data
Expand chromosome coordinates from the zero-based, half-open format.

```julia
using Bedgraph

const chrom_starts = [49302000, 49302300, 49302600, 49302900, 49303200, 49303500, 49303800, 49304100, 49304400]
const chrom_ends = [49302300, 49302600, 49302900, 49303200, 49303500, 49303800, 49304100, 49304400, 49304700]
const data_values = [-1.0, -0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75, 1.00]

(n, expanded_data_values) = Bedgraph.expand(chrom_starts, chrom_ends, data_values)
```

```julia

using Bedgraph

const records = [Record("chr19", 49302000, 49302300, -1.0), Record("chr19", 49302300, 49302600, -1.75)]

n, expanded_data_values = Bedgraph.expand(records)
```
