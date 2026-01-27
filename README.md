# Clusty

[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/clusty.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/clusty)
[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/clusty/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/clusty/releases)
[![Build and tests](../../workflows/Build%20and%20tests/badge.svg)](../../actions/workflows/main.yml)
[![License](https://anaconda.org/bioconda/famsa/badges/license.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

![x86-64](https://img.shields.io/static/v1?label=%E2%80%8B&message=x86-64&color=yellow&logo=PCGamingWiki&logoColor=white)
![ARM](https://img.shields.io/static/v1?label=%E2%80%8B&message=ARM&color=yellow&logo=Raspberry%20Pi&logoColor=white)
![Apple M](https://img.shields.io/static/v1?label=%E2%80%8B&message=Apple%20M&color=yellow&logo=Apple&logoColor=white)
![Windows](https://img.shields.io/badge/%E2%80%8B-Windows-00A98F?logo=windows)
![Linux](https://img.shields.io/static/v1?label=%E2%80%8B&message=Linux&color=00A98F&logo=linux&logoColor=white)
![macOS](https://img.shields.io/badge/%E2%80%8B-macOS-00A98F?logo=apple)


Clusty is a tool for large-scale clustering. By using sparse distance matrices it allows clustering data sets with millions of objects. 

## Major changes
- v1.3.0: Updated the Leiden graph construction, which affects cluster granularity. The previous bidirectional approach (IMG/VR-based; Camargo et al., 2023), where each connection contributed twice to the graph construction, was replaced with a single edge per connection. As a result, clustering at the same resolution value is finer than in earlier versions. To reproduce results from versions prior to v1.3.0, divide the value of `--leiden-resolution` by 2. 
- v1.2.0: Different ordering of assignments in the output file (clusters decreasingly by size, elements within clusters decreasingly by representativeness). Improved I/O performance.
- v1.1.0: Large memory optimizations. 

## Quick start

```bash
git clone --recurse-submodules https://github.com/refresh-bio/clusty
cd clusty
gmake -j

cd ./test

# Run single linkage clustering on the pairwise similarities stored in ictv.ani file, output cluster identifiers
# (two first columns are used as sequence identifiers, third one is assumed to store similarities).
../bin/clusty --algo single --similarity ictv.ani ictv.single

# Run uclust clustering accepting pairwise connectios with values greater or equal 0.70 in the ani column, output cluster representatives.
../bin/clusty --algo uclust --similarity --min ani 0.70 --out-representatives ictv.ani ictv.uclust.70

# Run CD-HIT clustering accepting pairwise connectios with values greater or equal 0.95 in the ani column, output cluster identifiers
# (use id2 and id2 columns as object identifiers and ani column as the similarity).
../bin/clusty --algo cd-hit --similarity --min ani 0.95 --id-cols id2 id1 --distance-col ani vir61.ani vir61.single.95

# Run complete linkage clustering, consider all objects from ictv.list file (including those without pairwise connections).
../bin/clusty --algo complete --objects-file ictv.list --similarity ictv.ani ictv.complete

```


## Installation

Clusty comes with a set of precompiled binaries for Windows, Linux, and macOS. They can be found under [Releases](./releases) tab.
The software is also available on [Bioconda](https://anaconda.org/bioconda/clusty):
```
conda install -c bioconda clusty
```
For detailed instructions how to set up Bioconda, please refer to the [Bioconda manual](https://bioconda.github.io/user/install.html#install-conda). 

The package can be built from the sources distributed as:
* Visual Studio 2022 solution for Windows,
* GNU Make project for Linux and macOS (gmake 4.3 and gcc/g++ 10 or newer required)

To compile Clusty under Linux/macOS please run:
```
gmake -j 
```

### Leiden algorithm

Clusty provides igraph's implementation of the Leiden algorithm. Precompiled binaries as well as bioconda distributions include Leiden algorithm. However, as igraph requires several external dependencies (CMake 3.18, Flex, Bison), it is by default not linked to the Clusty software. To install dependencies under Debian/Ubuntu linux use the following command:
```
sudo apt-get install cmake flex bison
```
Then, one needs to build the package with an additional option enabled:
```
gmake -j LEIDEN=true
```

Under Windows, Clusty is by default linked against igraph and it requires CMake as the only system dependency. After installing it (https://cmake.org) a user can run `build_igraph.bat` batch script which downloads Flex and Bison binaries to the appropriate locations and then builds igraph. After that it is possible to build Clusty using Visual Studio (the solution is located in `./src/clusty.sln`).


## Usage
`clusty [options] <distances> <assignments>`

Parameters:

* `<distances>` - input TSV/CSV table with pairwise distances
* `<assignments>` - output TSV/CSV file with cluster assignments

Options:

* `--objects-file <string>` - optional TSV/CSV file with object names in the first column sorted decreasingly w.r.t. representativness
* `--algo  <single | complete | uclust | set-cover | cd-hit | leiden>` - clustering algorithm:
  * `single`    - single linkage (connected component, MMseqs mode 1)
  * `complete`  - complete linkage
  * `uclust`    - UCLUST
  * `set-cover` - greedy set cover (MMseqs mode 0)
  * `cd-hit`    - CD-HIT (greedy incremental; MMseqs mode 2)
  * `leiden`    - Leiden algorithm
  
* `--id-cols <column-name1> <column-name2>` - names of columns with object identifiers (default: two first columns)
* `--distance-col <column-name>` - name of the column with pairwise distances (or similarities; default: third column)
* `--similarity` - use similarity instead of distance (has to be in [0,1] interval; default: false)
* `--percent-similarity` - use percent similarity (has to be in [0,100] interval; overrides `--similarity` flag, default: false)
* `--min <column-name> <real-threshold>` - accept pairwise connections with values greater or equal a given threshold in a specified column
* `--max <column-name> <real-threshold>` - accept only pairwise connections with values lower or equal a given threshold in a specified column
* `--numeric-ids` - use when sequences in the distances file are represented by numbers (can be mapped to string ids by the object file)
* `--out-representatives` - output representative objects for each cluster instead of cluster numerical identifiers
* `--out-csv` -- output a CSV table instead of a default TSV
* `-t` - number of threads (default: 4) 

Leiden algorithm options:

* `--leiden-resolution` - *resolution* parameter controlling clustering granularity in the Leiden algorithm (default: 0.7). Starting with v1.3.0, the underlying graph uses single edges instead of bidirectional connections. As a result, clustering at the same resolution value is coarser than in earlier versions. To reproduce results from versions prior to v1.3.0, use half of the resolution value previously applied.
* `--leiden-beta` - *beta* parameter for Leiden algorithm (default: 0.01)
* `--leiden-iterations` - number of interations for Leiden algorithm (default: 2)

## Examples

### Basic use case
The minimum input requirement is a TSV/CSV table with pairwise distances between objects (or similarities, if `--use-similarity` flag is used). By default, identifiers are assumed to be in the two first columns while distances are expected in the third one. Lack of a distance for a given pair of objects is translated to infinite distance. The example input table is given below:
```
name1	name2	ani
xxx	xx	0.93
aaa	aa	0.94
aaa	a	0.92
xx	x	0.94
bb	b	0.71
aa	a	0.89
b	bb	0.99
``` 

If the table is organized differenty, one can specify columns with ids and distances by using `--id-cols` and `--distance-col` parameters. The algorithm produces as a result a TSV table (or CSV if `--out-csv` flag is specified) with object identifiers followed by 0-based numerical identifiers of clusters. The clusters are ordered decreasingly by size with objects within a cluster outputed with decreasing representativeness (by default: in a lexigraphical order).
```
object	cluster
x	0
xx	0
xxx	0
a	1
aa	1
aaa	1
b	2
bb	2
```

### Cluster representatives
Instead of numerical identifiers, the package can output a representative object for every cluster using `--out-representatives` flag. The representative is the first object on the list that was assigned to the cluster:
```
object	cluster
x	x
xx	x
xxx	x
a	a
aa	a
aaa	a
b	b
bb	b
```
### Objects file
An optional TSV/CSV file specified by `--objects-file` parameter contains a complete list of analyzed objects. This file can be useful when distance table does not contain all the object (e.g., due to some filtering). Additionally, this file overrides representativeness of object in the ouput table. Thus, it affects ordering of objects within a cluster and their representative when `--out-representatives` flag is used. For instance, if one specifies the following objects file:
```
object
aaa
aa
a
bb
b
c
d
e
f
g
xxx
xx
x
```

the following output will be produced:
```
object	cluster
aaa	0
aa	0
a	0
xxx	1
xx	1
x	1
bb	2
b	2
c	3
d	4
e	5
f	6
g	7
```

or with `--out-representatives` flag:

```
object	cluster
aaa	aaa
aa	aaa
a	aaa
xxx	xxx
xx	xxx
x	xxx
bb	bb
b	bb
c	c
d	d
e	e
f	f
g	g

```

### Numerical identifiers

Clusty also supports numerical identifiers (`--numeric-ids` flag), e.g., as in table below. Note: the identifiers are assumed to be ordinal numbers (they do not have to be continous, though - gaps are allowed) rather than numerical hashes. This is because the software allocates memory area proportional to the size of the largest identifier.  
```
id1	id2	ani
10	11	0.93
0	1	0.94
0	2	0.92
11	12	0.94
3	4	0.71
1	2	0.89
4	3	0.99
5	6	0.33
```

The clustering produces the following output (lowest object id translates to the highest representativeness):
```
object	cluster
10	0
11	0
12	0
0	1
1	1
2	1
3	2
4	2
```

An external object file can also be used in this mode. It assumes, that the numerical identifier in the distance table translates to the 0-based index in the object file. The previously mentioned object file combined with `--out-representatives` flag produces the same table as the one generated without `--numeric-ids` switch:

```
object	cluster
aaa	aaa
aa	aaa
a	aaa
xxx	xxx
xx	xxx
x	xxx
bb	bb
b	bb
c	c
d	d
e	e
f	f
g	g
```

## Algorithms

In the following section one can find detailed information on clustering algorithms in Clusty, with *n* representing the number of objects (vertices) and *e* the number of distances (edges) in the data set (graph).


![clustering-steps](https://github.com/user-attachments/assets/08c1d9ee-94eb-49c2-9999-1c4850aa2cc9)


## Citation
Zielezinski A, Gudyś A, Barylski J, Siminski K, Rozwalak P, Dutilh BE, Deorowicz S. Ultrafast and accurate sequence alignment and clustering of viral genomes. Nat Methods. [https://doi.org/10.1038/s41592-025-02701-7](https://doi.org/10.1038/s41592-025-02701-7)

###
