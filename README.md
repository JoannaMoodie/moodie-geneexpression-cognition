# moodie-geneexpression-cognition
This repository contains the R scripts that I wrote in support of the preprint "General and specific patterns of cortical gene expression as spatial correlates of complex cognitive functioning". Find the preprint [here](https://www.biorxiv.org/content/10.1101/2023.03.16.532915v1). These scripts were run in R version 4.0.2. 

Please get in touch with me at jmoodie@ed.ac.uk if you have any quesitons.

## Brain regional morphometric profiles of _g_ 
Data sources: 
- These profiles are calculated with data from 3 cohorts, from which it is possible to request data from: [the UK Biobank](http://www.ukbiobank.ac.uk/register-apply/),  [the STratifying Resilience and Depression Longitudinally (STRADL) study](https://www.research.ed.ac.uk/en/datasets/stratifying-resilience-and-depression-longitudinally-stradl-a-dep) and the [Lothian Birth Cohort 1936](https://www.ed.ac.uk/lothian-birth-cohorts/data-access-collaboration).

## Brain regional profiles of gene expression
Data sources: 
- We used French and Paus' Deskian-Killiany cortical regional summary of the Allen Human Brain Atlas data [available here](https://figshare.com/articles/dataset/A_FreeSurfer_view_of_the_cortical_transcriptome_generated_from_the_Allen_Human_Brain_Atlas/1439749). 

### Validation of PCA components
Data sources: 
- The validation based on different pipelines for summarising the Allen Human Brain Atlas data uses data obtained with [Markello et al's scripts](https://github.com/netneurolab/markello_transcriptome) and the [abagen toolbox](https://github.com/rmarkello/abagen).

##  Correlations between _g_ and gene expression profiles
The cortical regional profiles (for gene expression and _g_) calculated in the previous scripts are correlated. The profiles of both general components (Component 1 and 2 from the PCA) and individual genes after controlling for these components are correlated with _g_-morphometry profiles. 

