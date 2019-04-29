#' @title GO_SS class object definition.
#' @description This class is invoked by \code{\link{build_GO_SS}} method in order to store \code{\link{enrich_GO_terms-class}} object, Information Content (IC),
#' and GO terms or groups distances objects based on semantic similarity.
#' @family GO_semantic_similarity
#' @slot db should be "Bioconductor", "EntrezGene", "Ensembl", or "Uniprot" ressource name.
#' @slot stamp date of stamp or annotation version if available in \code{character}.
#' @slot organism target species GO anotation in \code{character}.
#' @slot ont used ontology with "MF", "BP", or "CC".
#' @slot topGO  \code{list} with topGO objects summary informations.
#' @slot IC Information Content (IC)
#' @slot enrich_GOs \code{\link{merge_enrich_terms}} output object (\code{\link{enrich_GO_terms-class}} object).
#' @slot terms_dist \code{list} of GO terms or groups distances objects based on semantic similarity.
#' @include enrich_GO_terms.R
setClass(
    "GO_SS",
    slots=c(
        db="character",
        stamp ="character",
        organism="character",
        ont="character",
        topGO="list",
        IC ="numeric",
        enrich_GOs="enrich_GO_terms",
        terms_dist="list"
    )
)

#' @aliases GO_SS
#' @importFrom data.table melt.data.table .SD
setMethod(
    "show",
    signature="GO_SS",
    function(object){

        # Extract table
        Data<-slot(slot(object,"enrich_GOs"),"data")

        # Extract pvalues
        Data<-Data[,grep("\\.pvalue",names(Data)),with=FALSE]

        # count significant pvalues by condition
        Data<-Data[,lapply(.SD,function(x){sum(x<0.01)}),.SDcols=seq_len(ncol(Data))]

        # melt the table
        Data<-melt.data.table(
            Data,
            measure.vars=names(Data),
            variable.name = "conditions",
            value.name = "significant GO terms number"
        )

        # remove .pvalue in conditions column
        Data[,"conditions":=gsub("\\.pvalue","",Data$conditions)]

        # get topGO information
        topGO<-slot(object,"topGO")

        # format topGO information
        topGO<-vapply(names(topGO),function(x){

            # topGO subset
            Data<-topGO[[x]]

            # summery by element
            elems<-paste(
                vapply(names(Data),function(y){

                    # summery by element
                    paste(
                        paste("   ",y),
                        "\n       ",
                        paste(
                            paste(
                                names(Data[[y]]),
                                sub("\n.+$","",unlist(Data[[y]])),
                                sep=": "
                            ),
                            collapse="\n        "
                        )
                    )
                },""),
                collapse="\n  "
            )

            # summery by element
            paste(paste(x,elems,sep ="\n  "),"\n ")
        },"")

        # cat some text
        cat(
            "- object class: GO_SS",
            "\n- database: ",
            slot(object,"db"),
            "\n- stamp/version: ",
            slot(object,"stamp"),
            "\n- organism id: ",
            slot(object,"organism"),
            "\n- ontology: ",
            slot(object,"ont"),
            "\n- input:\n        ",
            paste(
                paste(
                    names(
                        slot(slot(object,"enrich_GOs"),"input")
                    ),
                    vapply(slot(slot(object,"enrich_GOs"),"input"),function(x){
                        paste(x,collapse=", ")
                    },""),
                    sep=": "
                ),
                collapse="\n        "
            ),
            "\n- topGO summary:\n ",
            topGO,
            "\n- enrich GOs data.table: ",
            nrow(
                slot(slot(object,"enrich_GOs"),"data")
            ),
            " GO terms of ",
            nrow(Data),
            " conditions.",
            paste(
                "\n       ",
                Data$conditions,
                ":",
                 Data$`significant GO terms number`,
                "terms"
            ),
            if(length(slot(object,"terms_dist"))>0){
                paste(
                    "\n- terms distances: ",
                    paste(
                        names(
                            slot(object,"terms_dist")
                        ),
                        collapse=", "
                    )
                )
            },
            sep=""
        )
    }
)
