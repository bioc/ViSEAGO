#' @title gene2GO class object definition.
#' @description This class is invoked by \code{\link{annotate}} method in order to store GO annotations for each category (MF, BP, CC).
#' @family GO_terms
#' @slot db database source in \code{character}.
#' @slot stamp date of stamp  in \code{character}.
#' @slot organism target species GO annotation  in \code{character}.
#' @slot MF a \code{list} containing  GO terms for Molecular Function (MF) category for each gene element.
#' @slot BP a \code{list} containing  GO terms for Biological Process (BP) category for each gene element.
#' @slot CC a \code{list} containing  GO terms for Cellular Component (CC) category for each gene element.
setClass(
    "gene2GO",
    slots=c(
        db="character",
        stamp="character",
        organism="character",
        MF="list",
        BP="list",
        CC="list"
    )
)

#' @aliases gene2GO
setMethod(
    "show",
    signature="gene2GO",
    function(object){

        # cat some text
        cat("- object class: gene2GO",
            "\n- database: ",slot(object,"db"),
            "\n- stamp/version: ",slot(object,"stamp"),
            "\n- organism id: ",slot(object,"organism"),
            "\n\nGO annotations:",
            "\n- Molecular Function (MF): ",length(slot(object,"MF"))," annotated genes with ",
            sum(lengths(slot(object,"MF")))," terms (",length(unique(unlist(slot(object,"MF")))),
            " unique terms)",
            "\n- Biological Process (BP): ",length(slot(object,"BP"))," annotated genes with ",
            sum(lengths(slot(object,"BP")))," terms (",length(unique(unlist(slot(object,"BP")))),
            " unique terms)",
            "\n- Cellular Component (CC): ",length(slot(object,"CC"))," annotated genes with ",
            sum(lengths(slot(object,"CC")))," terms (",length(unique(unlist(slot(object,"CC")))),
            " unique terms)",
            sep=""
        )
    }
)
