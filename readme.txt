Scripts for reproducing the analysis, figures and tables of the manuscript "Are
regular environmental spatial patterns periodic?" by Kastner et. al.  submitted
to Ecological Letters on 10 July 2023. This is the substantially revised and
extended first part of the unpublished manuscript with working title
"Irregularity of self-organized vegetation patterns: Theory and empirical
evidence" submitted to PNAS on 19 January 2022.

Requirement:
- GNU/Linux (recommended operating system)
- Python 3
- Matlab R2022b (older versions probably work, toolbox dependencies listed in
                 in pattern_periodicity_test_metadata.m)
- dependencies from https://github.com/karlkastner

The script "pattern_periodicity_test_batch" reproduces the whole workflow.
It first fetches dependencies from GitHub, and subsequently executes scripts 
for fetching and analyzing the patterns, followed by plotting the results.

Processing the whole database (~10000) patterns takes approximately 24-48 h.
For quick testing, set the parameter skip in pattern_periodicity_test_metadata.m
to a large integer such as 10, 100 or 1000. In this case only every skip pattern
will be processed.

The GPL v. 3 extends over all files in the repository, i.e. computer scripts,
documentation and geospatial data files, except for the sample satellite images
in the input folder, which are property of Google / Maxar Technologies (2023).  

