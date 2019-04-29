#' @title enrich_GO_terms class object definition.
#' @description This class is invoked by \code{\link{merge_enrich_terms}} method in order to store the merged \code{data.table} and associated metadata.
#' @family enrich_GO_terms
#' @slot same_genes_background logical.
#' @slot input a list containing named elements. Each element must contain the name of \code{\link[topGO]{topGOdata-class}}
#' object created by \code{\link{create_topGOdata}} method and the associated  \code{\link[topGO]{topGOresult-class}}
#' object(s) to combinate (see examples in \code{\link{merge_enrich_terms}}).
#' @slot ont ontology used "MF", "BP", or "CC".
#' @slot topGO  a \code{list} with topGO objects summary informations.
#' @slot data a merged \code{data.table}  of enriched GO terms (p<0.01) in at least once with GO descriptions and statistical values.
setClass(
    "enrich_GO_terms",
    slots=c(
        same_genes_background="logical",
        input="list",
        ont="character",
        topGO="list",
        data="data.table"
    )
)

#' @aliases enrich_GO_terms
#' @importFrom data.table .SD melt.data.table :=
setMethod(
    "show",
    "enrich_GO_terms",
    function(object){

        # Extract table
        Data<-slot(object,"data")

        # Extract pvalues
        Data<-Data[,grep("\\.pvalue",names(Data)),with=FALSE]

        # count significant pvalues by condition
        Data<-Data[,lapply(.SD,function(x){sum(x<0.01,na.rm = TRUE)}),.SDcols=seq_len(ncol(Data))]

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
        cat("- object class: enrich_GO_terms",
            "\n- ontology: ",slot(object,"ont"),
            "\n- input:\n        ", paste(paste(names(object@input),
            vapply(slot(object,"input"),function(x){paste(x,collapse=", ")},""),sep=": "),collapse="\n        "),
            "\n- topGO summary:\n ", topGO,
            "\n- enrich GOs data.table (p<0.01 in at least one list): ",nrow(slot(object,"data"))," GO terms of ",nrow(Data)," conditions.",
            paste("\n       ",Data$conditions,":",Data$`significant GO terms number`,"terms"),sep=""
        )
    }
)
