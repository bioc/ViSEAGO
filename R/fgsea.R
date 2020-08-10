#' @title fgsea class object definition.
#' @description This class is invoked by \code{\link{runfgsea}} method in order to store results.
#' @slot description a \code{character} string with database source, date of stamp, and target species GO annotation.
#' @slot method fgsea method used.
#' @slot params a \code{list} containing used input parameters for perform \code{\link[fgsea]{fgseaSimple}} or \code{\link[fgsea]{fgseaMultilevel}}.
#' @slot input a \code{list} containing input values.
#' @slot data a \code{list} containing \code{data.table} fgsea procedure output.
setClass(
    "fgsea",
    slots=c(
        description="character",
        ontology="character",
        method="character",
        params="list",
        input="list",
        data="list"
    )
)

#' @aliases fgsea
setMethod(
    "show",
    signature="fgsea",
    function(object){

        # extract database informations
        db=unlist(
            strsplit(
                slot(object,"description"),
                " "
            )
        )

        # get parameters informations
        params<-paste(
            "\n    ",
            names(
                slot(object,"params")
            ),
            " : ",
            unlist(
                as.character(slot(object,"params"))
            ),
            sep=""
        )
        
        # get fgsea results
        data<-slot(object,"data")[[1]]

        # cat some text
        cat("- object class: fgsea",
            "\n- database: ",db[1],
            "\n- stamp/version: ",db[3],
            "\n- organism id: ",db[2],
            "\n- ontology: ",slot(object,"ontology"),
            "\n- method: ",slot(object,"method"),
            "\n- fgseaMultilevel parameters: ",params,
            "\n- fgseaMultilevel results:",
            "\n    ",nrow(data), " tested GO terms with enrichment pvalue from ",min(data$pval)," to ",max(data$pval),
            sep=""
        )
    }
)
