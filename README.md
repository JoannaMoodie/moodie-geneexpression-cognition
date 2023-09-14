# moodie-geneexpression-cognition
This repository contains the R scripts that I wrote in support of the preprint "General and specific patterns of cortical gene expression as spatial correlates of complex cognitive functioning" [available here](https://www.biorxiv.org/content/10.1101/2023.03.16.532915v1). These scripts were run in R version 4.0.2. 

Please get in touch with me at jmoodie@ed.ac.uk if you have any quesitons.

## scripts
### Brain regional morphometric profiles of _g_ 
Script: brainregion_g_morphometry.R

Data sources: 
- These profiles are calculated with data from 3 cohorts, from which it is possible to request data from: [the UK Biobank](http://www.ukbiobank.ac.uk/register-apply/),  [the STratifying Resilience and Depression Longitudinally (STRADL) study](https://www.research.ed.ac.uk/en/datasets/stratifying-resilience-and-depression-longitudinally-stradl-a-dep) and the [Lothian Birth Cohort 1936](https://www.ed.ac.uk/lothian-birth-cohorts/data-access-collaboration).

### Brain regional profiles of gene expression
Script: brainregion_geneexpression.R

Data sources: 
- We used [French and Paus' 2015](https://figshare.com/articles/dataset/A_FreeSurfer_view_of_the_cortical_transcriptome_generated_from_the_Allen_Human_Brain_Atlas/1439749) Deskian-Killiany cortical regional summary of the Allen Human Brain Atlas data. 
- The validation of components based on different pipelines for summarising the Allen Human Brain Atlas data uses data obtained with [Markello et al's scripts](https://github.com/netneurolab/markello_transcriptome) and the [abagen toolbox](https://github.com/rmarkello/abagen).
- The cell types were categorized according to [Zeisel et al.'s (2015) paper](https://pubmed.ncbi.nlm.nih.gov/25700174/).

###  Correlations between _g_ and gene expression profiles
Script: brainregion_correlations_g_geneexpression.R

Data source:
- The cortical regional profiles (for gene expression and _g_) calculated in the previous scripts are correlated.

## data
- UKB_neuroexclusions.csv is a list of UKB IDs that we excluded on the basis of neurological conditions (see manuscript for details).
- component_loadings.csv is a list of component loadings after varimax rotation for Component 1 and Component 2 (for 8235 genes).
- regional_profiles.xlsx contains the the regional profiles of the two gene expression components and the regional meta-analysis results for g-morphometry associations. 
