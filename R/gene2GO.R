#' @title gene2GO class object definition.
#' @description This class is invoked by \code{\link{annotate}} method in order to store GO annotations for each category (MF, BP, CC).
#' @importFrom methods setClass slot
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
setMethod("show",signature="gene2GO",function(object){

  ###################
  # cat some text
  base::cat("- object class: gene2GO",
    "\n- database: ",methods::slot(object,"db"),
    "\n- stamp/version: ",methods::slot(object,"stamp"),
    "\n- organism id: ",methods::slot(object,"organism"),
    "\n\nGO annotations:",
    "\n- Molecular Function (MF): ",base::length(methods::slot(object,"MF"))," annotated genes with ",
    base::sum(base::lengths(methods::slot(object,"MF")))," terms (",base::length(base::unique(base::unlist(methods::slot(object,"MF")))),
    " unique terms)",
    "\n- Biological Process (BP): ",base::length(methods::slot(object,"BP"))," annotated genes with ",
    base::sum(base::lengths(methods::slot(object,"BP")))," terms (",base::length(base::unique(base::unlist(methods::slot(object,"BP")))),
    " unique terms)",
    "\n- Cellular Component (CC): ",base::length(methods::slot(object,"CC"))," annotated genes with ",
    base::sum(base::lengths(methods::slot(object,"CC")))," terms (",base::length(base::unique(base::unlist(methods::slot(object,"CC")))),
    " unique terms)",sep="")
})
