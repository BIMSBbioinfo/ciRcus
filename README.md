# ciRcus: an R package for circRNA manipulation, annotation and analysis

ciRcus is a collection of functions for everyday munging of circRNA data.
In its current, preliminary version, it can take lists of putative splice junctions generated using find_circ.py ([Memczak et al. 2013](http://www.nature.com/nature/journal/v495/n7441/full/nature11928.html), [circBase](http://www.circbase.org)) as input, and perform following annotation steps:

* quality filtering
* suggesting a host gene candidate was spliced from
* calculating circular-to-linear ratio
* describe gene features circRNA candidate was spliced from
* report if circRNA's spliced junctions are already annotated as linear exon-intron junctions
* generate read count histogram and gene feature pie-charts

# Installation

### Install via install_github()
```R
#' Install dependecies
install.packages( c("data.table", "DBI", "hash", "ggplot2", "RMySQL", "devtools"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges","GenomicFeatures", "IRanges", "biomaRt", "AnnotationHub"))

#' install the package
library(devtools)
install_github("BIMSBbioinfo/ciRcus", build_vignettes=FALSE)


```

# Using the package
### Build TxDb object with genomic features and save locally
Load genomic features from Ensembl and build a database for later (re)use. Currently supported assemblies are hg19, hg38, mm10, rn5 and dm6. This needs to be done only once per assembly. 
```R
gtf2sqlite( assembly = "hg19",
            db.file  = system.file("~/human_hg19_ens75_txdb.sqlite",
                                   package="ciRcus")
```
### Extract features from the database
List of features returned by `loadAnnotation()` will be used to annotate circRNAs. Saving it as a separate object is a good practice once we start analyzing multiple circRNA libraries.
```R
annot.list <- loadAnnotation(system.file("~/human_hg19_ens75_txdb.sqlite",
                                         package="ciRcus"))
```
### Load and annotate circRNAs
```R
cdata <- data.frame(sample=c("FC1", "FC2", "H1", "H2", "L1", "L2"),
                    filename=list.files(system.file('extdata', package='ciRcus'),                                                                         pattern='sites.bed',
                                        full.names=TRUE)[1:6])
circs.se <- summarizeCircs(colData=cdata, wobble=1)
circs.se <- annotateCircs(circs.se, annot.list=annot.list)
circs.dt <- resTable(circs.se)
circs.dt
```
### Plot data
```R
histogram(circs.se, 0.5)
annotPie(circs.se, 0.02)
```
