#' @title genomic_ressource class object definition.
#' @description This class stores the annotations and associated metadata obtained by \code{\link{Bioconductor2GO}},
#' \code{\link{EntrezGene2GO}}, \code{\link{Ensembl2GO}}, or \code{\link{Uniprot2GO}} methods.
#' @importFrom methods setClass slot
#' @family genomic_ressource
#' @slot db name of database used (Bioconductor, EntrezGene, Ensembl, or Uniprot).
#' @slot stamp date of stamp (for Bioconductor, EntrezGene, and Uniprot), or annotation version for Ensembl database.
#' @slot data GO annotations from \code{\link{EntrezGene2GO}} method.
#' @slot organisms informations about species/datasets availables.
#' @slot mart Ensembl mart from \code{\link{Ensembl2GO}} method.
setClass("genomic_ressource",
         slots=c(
            db="character",
            stamp="character",
            data="data.table",
            organisms="data.table",
            mart="list"
        )
)
#' @importFrom methods setMethod
setMethod("show", signature="genomic_ressource",function(object){

  ###################
  # cat some text
  base::cat("- object class: genomic_ressource",
    "\n- database: ",methods::slot(object,"db"),
    "\n- stamp/version: ",methods::slot(object,"stamp"),
    "\n- available organisms: ",base::nrow(methods::slot(object,"organisms")),
  sep="")
})
