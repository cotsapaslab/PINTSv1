 
#### Description of the package

The PINTS (Protein Interaction Network Tissue Search) package is designed to provide a framework to identify a group of 
interacting genes under mutational constraint and equivalently, a group of interacting genes associated with diseases/traits 
(e.g., trait-associated subnetwork from genome wide association study) and to infer the biological context represented by the
group of genes, e.g., likely tissue of action. These can help to gain an insight about the molecular mechanism associated with 
mutational constraint and complex diseases/traits. 


#### Installation

You can install the PINTS package from github.

        install.packages('devtools')
        library(devtools)
        
        install_github("PINTS", username="PINTS-CotsapasLab")
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


#### Acknowledgement

We thank all the members of Cotsapas Lab for the develpoment of the project and also thank the whole research community for 
making the useful datasets and tools available. 
