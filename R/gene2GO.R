#' @title gene2GO class object definition.
#' @description This class is invoked by \code{\link{annotate}} method in order to store GO annotations for each category (MF, BP, CC).
#' @importFrom methods setClass
#' @family GO_terms
#' @slot db database source in \code{character}.
#' @slot stamp date of stamp  in \code{character}.
#' @slot organism target species GO annotation  in \code{character}.
#' @slot MF a \code{list} containing  GO terms for Molecular Function (MF) category for each gene element.
#' @slot BP a \code{list} containing  GO terms for Biological Process (BP) category for each gene element.
#' @slot CC a \code{list} containing  GO terms for Cellular Component (CC) category for each gene element.
setClass("gene2GO",
         slots=c(
           db="character",
           stamp="character",
           organism="character",
           MF="list",
           BP="list",
           CC="list"
         )
)
#' @importFrom methods setMethod
setMethod("show", "gene2GO",function(object) {

  ###################
  # cat some text
  base::cat("- object class: gene2GO",
    "\n- database: ",object@db,
    "\n- stamp/version: ",object@stamp,
    "\n- organism id: ",object@organism,
    "\n\nGO annotations:",
    "\n- Molecular Function (MF): ",base::length(object@MF)," annotated genes with ",
    base::sum(base::lengths(object@MF))," terms (",base::length(base::unique(base::unlist(object@MF))),
    " unique terms)",
    "\n- Biological Process (BP): ",base::length(object@BP)," annotated genes with ",
    base::sum(base::lengths(object@BP))," terms (",base::length(base::unique(base::unlist(object@BP))),
    " unique terms)",
    "\n- Cellular Component (CC): ",base::length(object@CC)," annotated genes with ",
    base::sum(base::lengths(object@CC))," terms (",base::length(base::unique(base::unlist(object@CC))),
    " unique terms)",sep="")
})
