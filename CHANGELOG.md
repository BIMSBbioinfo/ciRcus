ciRcus 0.1.7
------------
IMPROVEMENTS AND BUG FIXES

* 18-column find_circ output (wherever it comes from) supported

ciRcus 0.1.6
------------
NEW FEATURES

* new bedTracks()/writeBedTracks() functions generating BED6-style
  GRangesList/writing BED files for circRNA browser tracks
* added D. melanogaster (dm6) and R. norvegicus (rn6) support

IMPROVEMENTS AND BUG FIXES

* using "coarse" pruning.mode in keepStandardChromosomes when building annotation
* introduced a number of bugfixes (and probably new bugs)

ciRcus 0.1.5
------------
NEW FEATURES

* loadAnnotation() converted to S4, can digest both TxDb objects and their sqlite dumps

IMPROVEMENTS AND BUG FIXES

* circLinRatio bug that prevented the ratio from being calculated in 0.1.4 fixed
* many under-the-hood improvements (travis, lintr, codecov...) users will indirectly benefit from

ciRcus 0.1.4
------------

NEW FEATURES

* CIRI2 support
* WBcel235 genome assembly supported
* QC plot for total vs. unique reads

IMPROVEMENTS AND BUG FIXES

* sample order from colData is not broken by summarizeCircs
* system.file call fixed in README


ciRcus 0.1.3
------------
NEW FEATURES

* find_circ2 supported

IMPROVEMENTS AND BUG FIXES

* different coordinate indexing in input and annotation
  will be automatically detected and fixed
* lots of under-the-hood improvements

ciRcus 0.1.2
------------
NEW FEATURES

* find_circ.py v1.2 supported

ciRcus 0.1.1
--------------
NEW FEATURES

* input is loaded into a SummarizedExperiment object
* basic circBase interaction

  
IMPROVEMENTS AND BUG FIXES

* unit tests
