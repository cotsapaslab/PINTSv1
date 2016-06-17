 
#### Description of the package

The PINTS (Protein Interaction Network Tissue Search) package provides a framework to identify groups of interacting genes with high user--specified scores and perform tissue--specific expression enrichment for any significant groups identified. We have (used it to demonstrate)[http://dx.doi.org/10.1371/journal.pgen.1006121 "PLoS Genetics PINTS paper"] that genes under purifying selection in the human population are clustered, suggesting selection acts on entire biological mechanisms.


#### Installation

You can install the PINTS package from github.

        install.packages('devtools')
        library(devtools)
        
        install_github("CotsapasLab/PINTSv1")
        library(PINTS)

In addition to this, you would need to install all package dependencies. They are downloadble frm CRAN 
(http://cran.r-project.org) and Bioconductor (http://www.bioconductor.org) project website. 
The following are the necessary packages.

        #Installing from CRAN
        install.packages("igraph")

        #installing from Bioconductor.org
        source("http://bioconductor.org/biocLite.R")
        biocLite("BioNet")
        biocLite("RBGL")  
        biocLite("graph")  
        biocLite("Matching")

For the tissue-level context search in the step 3 of the PINTS workflow, you also need to have a UGM software running under 
MATLAB. The UGM software is available from http://www.di.ens.fr/~mschmidt/Software/UGM.html. We provide a set of Matlab scripts 
used for PINTS workflow separately (Download inst/extdata/UGM_2011_CotsapasAdd.zip and upzip the file under a working direcotry). 


#### Tutorials

A series of tutorials is intended to demonstrate the PINTS workflow to identify the mutational constraint gene subnetwork and 
to analyize the biological context (tissue specificity) of the identified disease-associated subnetwork.

You can check out the tutorials in a R session after the installation of the package.

        browseVignettes("PINTS")


#### Citation

Please cite (JinMyung Choi et al)[http://dx.doi.org/10.1371/journal.pgen.1006121], Network Analysis of Genome-Wide Selective Constraint Reveals a Gene Network Active in Early Fetal Brain Intolerant of Mutation, PLoS Genetics 12(6):e1006121 2016.
