# Code_AAV_integrations
_AAV integrations using RAAVioli in Ferrari, Cesana et al, Cell Stem Cell_

The code here released supported the generation of the results presented in the manuscript from Ferrari S., Jacob A., Cesana D., et al 2022.

## Pre-requisites
1. Install the pipeline [VISPA2](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1937-9) (Spinozzi G., et al 2017) from [this repository](https://github.com/giuliospinozzi/vispa2): https://github.com/giuliospinozzi/vispa2 
2. To run the next steps, you will need: Python 3, the Python libraries pandas and pysam.

## Instructions
Here we divided the code in the main steps and we structured the repository file system with the 3 main folders related the the steps.
Before running the code, you will need to create an hybrid genome (mixing your target genome, such as human reference hg38, with the rAAV genome), and indexing the hybrid genome with BWA.

### Step 1
Perform the alignment of the raw reads using VISPA2. The scratch code is reported in the folder step1. All quality procedures and barcode split/de-multiplexing are included in the code.

### Step 2
The Python code supports the identification of integration sites (IS) and vector-genome junctions.

### Step 3
Generation of the final matrix of IS with the corresponding quantification by number of reads (sequence count) and number of shearing sites (number of genomic fragments). The matrix is composed as follows: each row is an independent IS, labeld with the genomic annotations (such as genomic alignment start, corresponding to the IS by experimental design, and strand); in columns, the matrix shows the different source samples analyzed, containing the observed IS. Each cell contains the number of sequencing reads (coverage) for the observed IS or the number of different genomes (accounting the number of distinct genomic fragments).
