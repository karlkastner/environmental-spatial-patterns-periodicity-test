Source code and input data for reproducing the analysis of the manuscript:
"Unravelling the spatial structure of regular environmental spatial patterns",
submitted to Catena on 19 Feb 2024 (CATENA25913). Originally submitted with
working title "Are Environmental Spatial Patterns Periodic" to Catena
on 26 Sep 2023 (CATENA23872). Evolved from first part of the manuscript
"Irregularity of self-organized vegetation patterns: Theory and empirical
evidence", submitted to PNAS on 19 Jan 2022 (2022-00974).

Requirement:
- GNU/Linux (recommended operating system)
- Python 3
- Matlab R2023b (older versions probably work, toolbox dependencies listed in
                 in pattern_periodicity_test_metadata.m)
- dependencies from https://github.com/karlkastner

A snapshot of this repository which already includes the required library files
is available on: https://zenodo.org/records/13695204

The script "pattern_periodicity_test_batch" reproduces the whole workflow.
It first fetches dependencies from several GitHub repositories and subsequently
executes scripts:
- fetching and installing required library files,
- fetching the satellite images from the google maps tile server,
- analyzing the images,
- tabulating and plotting the results.

Processing the whole database (~10000) patterns takes several days.
For quick testing, set the parameter skip in pattern_periodicity_test_metadata.m
to a large integer such as 10, 100 or 1000. In this case only every skip pattern
will be processed.

The GPL v. 3 extends over all files in the repository, i.e. computer scripts,
documentation and geospatial data files, except for the sample satellite images
in the input folder, which are property of Google / Maxar Technologies (2023).  

