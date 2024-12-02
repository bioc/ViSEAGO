#' @title ViSEAGO package
#' @description Easier data mining of biological functions organized into clusters using Gene Ontology and semantic.
#' @details The main objective of ViSEAGO workflow is to carry out a data mining of biological functions and establish
#' links between genes involved in the study. We developed ViSEAGO in R to facilitate functional Gene Ontology (GO)
#' analysis of complex experimental design with multiple comparisons of interest.
#'
#' It allows to study large-scale datasets together and visualize GO profiles to capture biological knowledge.
#' The acronym stands for three major concepts of the analysis: Visualization, Semantic similarity and Enrichment Analysis of Gene Ontology
#' (\code{\link{pkgdiagram}}).
#'
#' It provides access to the last current GO annotations (\code{\link{annotate}}), which are retrieved from one of
#' NCBI EntrezGene (\code{\link{Bioconductor2GO}}, \code{\link{EntrezGene2GO}}),
#' Ensembl (\code{\link{Ensembl2GO}}) or Uniprot (\code{\link{Uniprot2GO}}) databases
#' for available species (\code{\link{available_organisms}}).
#'
#' ViSEAGO extends classical functional GO analysis (\code{\link{create_topGOdata}}) to focus on functional coherence
#' by aggregating closely related biological themes while studying multiple datasets at once (\code{\link{merge_enrich_terms}}).
#'
#' It provides both a synthetic and detailed view using interactive functionalities respecting the GO graph structure
#' (\code{\link{MDSplot}}, \code{\link{GOterms_heatmap}}, \code{\link{GOclusters_heatmap}}), and ensuring functional
#' coherence supplied by semantic similarity (\code{\link{build_GO_SS}}, \code{\link{compute_SS_distances}}).
#'
#' ViSEAGO has been successfully applied on several datasets from different species with a variety
#' of biological questions. Results can be easily shared between bioinformaticians and biologists, enhancing reporting capabilities while
#'  maintaining reproducibility.
#'
#' @docType _PACKAGE
#' @name ViSEAGO
NULL
