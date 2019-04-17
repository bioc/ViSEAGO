<img src="./inst/extdata/univ.png" align="right"/>
<img src="./inst/extdata/boa.png" align="right"/>
<img src="./inst/extdata/inra.png" align="right"/>

# ViSEAGO: Easier data mining of biological functions annotated with Gene Ontology and clusterized with semantic similarity

The main objective of ViSEAGO workflow is to carry out a data mining of biological functions and establish links between genes involved in the study. We developed ViSEAGO in R to facilitate functional Gene Ontology (GO) analysis of complex experimental design with multiple comparisons of interest. It allows to study large-scale datasets together and visualize GO profiles to capture biological knowledge. It stands for three major concepts: Visualization, Semantic similarity and Enrichment Analysis of Gene Ontology. It provides access to the last current GO annotations, which are retrieved from one of NCBI EntrezGene, Ensembl or Uniprot databases for available species. ViSEAGO extends classical functional GO analysis to focus on functional coherence by aggregating closely related biological themes while studying multiple datasets at once. It provides both a synthetic and detailed view using interactive functionalities respecting the GO graph structure and ensuring functional coherence supplied by semantic similarity. ViSEAGO has been successfully applied on several datasets from different species with a variety of biological questions. Results can be easily shared between bioinformaticians and biologists, enhancing reporting capabilities while maintaining reproducibility.

## Citation

to complete

## Installation

```r
###################
# load ViSEAGO package from Bioconductor
BiocManager::install("ViSEAGO")

###################
# load ViSEAGO package from gitlab
remotes::install_gitlab(
  "umr-boa/viseago",
  host = "forgemia.inra.fr",
  build_opts = c("--no-resave-data","--no-manual")
)
```

## Overview

![](./inst/extdata/figure2_BMCBioData.png)
