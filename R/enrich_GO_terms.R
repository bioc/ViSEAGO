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

        # pvalues threshlod according condition
        p<-vapply(slot(object,"topGO"),function(x){

            #  unlist
            x=unlist(x)

            # extract pvalue threshold
            as.numeric(
                sub(
                    "^.+<",
                    "",
                    x[grep("test_name",names(x))]
                )
            )
        },0)

        # count significant pvalues by condition
        Data<-lapply(seq_len(ncol(Data)),function(x){
            
            data.table(
                conditions=sub("\\.pvalue","",names(Data)[x]),
                `significant GO terms number`=sum(Data[,x,with=FALSE]<p[x],na.rm=TRUE)
            )
        })
        
        # bind results
        Data<-rbindlist(Data)

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
            "\n- enrich GOs (in at least one list): ",nrow(slot(object,"data"))," GO terms of ",nrow(Data)," conditions.",
            paste("\n       ",Data$conditions,":",Data$`significant GO terms number`,"terms"),sep=""
        )
    }
)
