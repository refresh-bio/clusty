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

## Quick start

```bash
git clone --recurse-submodules https://github.com/refresh-bio/clusty
cd clusty
make -j

cd ./test

# Run single linkage clustering on the pairwise similarities stored in ictv.ani file, output cluster identifiers
# (two first columns are used as sequence identifiers, third one is assumed to store similarities).
../clusty --algo single --similarity ictv.ani ictv.single

# Run uclust clustering accepting pairwise connectios with values greater or equal 0.70 in the ani column, output cluster representatives.
../clusty --algo uclust --similarity --min ani 0.70 --out-representatives ictv.ani ictv.uclust.70

# Run CD-HIT clustering accepting pairwise connectios with values greater or equal 0.95 in the ani column, output cluster identifiers
# (use id2 and id2 columns as object identifiers and ani column as the similarity).
../clusty --algo cd-hit --similarity --min ani 0.95 --id-cols id2 id1 --distance-col ani vir61.ani vir61.single.95

# Run complete linkage clustering, consider all objects from ictv.list file (including those without pairwise connections).
../clusty --algo complete --objects-file ictv.list --similarity ictv.ani ictv.complete

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
* MAKE project for Linux and macOS (g++-10 required).

To compile Clusty under Linux/macOS please run:
```
make -j 
```

### Leiden algorithm

Clusty provides igraph's implementation of the Leiden algorithm. Precompiled binaries as well as bioconda distributions include Leiden algorithm. However, as igraph requires several external dependencies (CMake 3.18, Flex, Bison), it is by default not linked to the Clusty software. To install dependencies under Debian/Ubuntu linux use the following command:
```
sudo apt-get install cmake flex bison
```
Then, one need to build the package with an additional option enabled:
```
make -j LEIDEN=true
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

* `--leiden-resolution` - resolution parameter for Leiden algorithm (default: 0.7)
* `--leiden-beta` - beta parameter for Leiden algorithm (default: 0.01)
* `--leiden-iterations` - number of interations for Leiden algorithm (default: 2)


The minimum input requirement is a CSV/TSV table with pairwise distances between objects (or similarities, if `--use-similarity` flag is used). By default, identifiers are assumed to be in the two first columns while distances are expected in the third one. Lack of a distance for a given pair of objects is translated to infinite distance. The example input table is given below:
```
id1,id2,distance
a,b,0.04
a,f,0.25
a,c,0.02
d,b,0.70
b,c,0.51
e,f,0.01
``` 

If the table is organized differenty, one can specify columns with ids and distances by using `--id-cols` and `--distance-col` parameters. The algorithm produces as a result a TSV table (or CSV if `--out-csv` flag is specified) with object identifiers followed by 0-based numerical identifiers of clusters. The objects appear in the order as in the distance table.
```
object	cluster
a	0
b	0
f	1
c	0
d	2
e	1
```

Alternatively, instead of numerical identifiers, package can output a representative object for every cluster using `--out-representatives` flag. The representative of a cluster is the first object on the list that was assigned to the cluster:
```
object	cluster
a	a
b	a
f	f
c	a
d	d
e	f
```

An optional TSV/CSV file specified by `--objects-file` parameter contains a complete list of analyzed objects. This file can be useful when distance table does not contain all the object (e.g., due to some filtering). Additionally, this file determines order in which objects appear in ouput table (thus, it also affects cluster representatives when using `--out-representatives` flag). For instance, if one specifies following objects file:
```
objects
b
x
a
c
e
d
f
```

the following output will be produced:
```
object	cluster
b	b
x	x
a	b
c	b
e	e
d	d
f	e
```
## Algorithms

In the following section one can find detailed information on clustering algorithms in Clusty, with *n* representing the number of objects (vertices) and *e* the number of distances (edges) in the data set (graph).

| Algorithm  | Details | Time complexity |
| ------------- | ------------- | ------------- |
| Single linkage | Hierarchical agglomerative clustering with a distance between groups defined as a distance between their closest members. Equivalent to finding all consistent subgraphs in a graph. Performed using breadth-first search. | *O*(*e*) |
| Complete linkage | Hierarchical agglomerative clustering with a distance between groups defined as a distance between their furthest members. Equivalent to finding a disjoint set of complete subgraphs covering the entire graph. Fast identification of clusters to merge is performed by storing distances in a heap. | *O*(*e* log*e*) |
| UCLUST | Greedy clustering with objects investigated in descending order w.r.t. representativeness. The first object becomes a centroid; the following objects are either (a) assigned to the closest of centroids they are connected with or (b) become new centroids if they are not connected to any of the existing ones. | *O*(*e*) |
| Greedy set cover | Greedy clustering with objects investigated in descending order w.r.t. the number of neighbors. Every unassigned object becomes a new centroid with all connected objects being assigned to it. Equivalent to MMseqs mode 0 clustering.  |  *O*(*n* log*n* + *e*) |
| CD-HIT | Greedy clustering with objects investigated in descending order w.r.t. the representativeness. Every unassigned object becomes a new centroid with all connected objects being assigned to it. Equivalent to MMseqs mode 2 clustering. | *O*(*e*) |
| Leiden | Iterative heuristic for finding communities in networks. It consists of three phases: (1) local moving of nodes, (2) refinement of the partition, and (3) aggregation of the network using the refined and the non-refined partitions. | unknown (algorithm provided by an external library) |

## Citation
Zielezinski A, Gudyś A, Barylski J, Siminski K, Rozwalak P, Dutilh BE, Deorowicz S. Ultrafast and accurate sequence alignment and clustering of viral genomes. bioRxiv [doi:10.1101/2024.06.27.601020].

###
