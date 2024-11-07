#pragma once

// *******************************************************************************************
// This file is a part of Clusty software distributed under GNU GPL 3 license.
// The homepage of the Clusty project is https://github.com/refresh-bio/Clusty
//
// Copyright(C) 2024-2024, A.Gudys, K.Siminski, S.Deorowicz
//
// *******************************************************************************************

#define VERSION		"1.1.4"
#define DATE		"2024-11-06"


/********* Version history *********

1.1.4 (2024-11-06)
* macOS compilation fix: SDK setting for igraph.

1.1.3 (2024-10-21)
* Added `--version` switch.

1.1.2 (2024-10-13)
* Precompiled binaries for macOS include Leiden algorithm.
* Fixed small bug with `--leiden-iterations` param being displayed in help as `--leiden-resolution`.

1.1.1 (2024-10-05)
Fixed bug with Leiden algorithm always running with default parameters.

1.1.0 (2024-09-30)
Memory optimizatons:
* Row identifier removed from dist_t structure.
* Object identifiers in complete linkage stored as 32-bit integers.
* For some clustering algorithms only edges are stored (without distances). 

1.0.3 (2024-09-12)
Fixes in building scripts.

1.0.2 (2024-08-29)
Fixed crash in reading very large input matrices. Reduced memory requirements.

1.0.1 (2024-08-19)
Fixed bug in parsing floating point numbers with more than 15 decimal places.

1.0.0 (2024-06-26)
First public release

*/
