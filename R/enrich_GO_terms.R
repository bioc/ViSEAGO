#' @title enrich_GO_terms class object definition.
#' @description This class is invoked by \code{\link{merge_enrich_terms}} method in order to store the merged \code{data.table} and associated metadata.
#' @family enrich_GO_terms
#' @slot same_genes_background logical.
#' object(s) to combinate (see examples in \code{\link{merge_enrich_terms}}).
#' @slot ont ontology used "MF", "BP", or "CC".
#' @slot method enrichment test used "topGO", or "fgsea".
#' @slot summary  a \code{list} with topGO or fgsea object(s) summary informations.
#' @slot data a merged \code{data.table}  of enriched GO terms (p<0.01) in at least once with GO descriptions and statistical values.
setClass(
    "enrich_GO_terms",
    slots=c(
        same_genes_background="logical",
        ont="character",
        method="character",
        cutoff="list",
        summary="list",
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
        p<-unlist(slot(object,"cutoff"))

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
        obj_summary<-slot(object,"summary")

        # summary according method
        obj_summary<-vapply(names(obj_summary),function(x){

            # subset
            Data<-obj_summary[[x]]

            # topGO summary
            if(slot(object,"method")=="topGO"){

                # topGO summary by element
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

                # summary by element
                return(paste(paste(x,elems,sep ="\n  "),"\n "))
            }

            # fgsea summary
            if(slot(object,"method")=="fgsea"){

                # fgsea summary by element
                elems<-paste(
                    paste(
                        "\n    ",
                        names(
                            Data
                        ),
                        " : ",
                        unlist(
                            Data
                        ),
                        sep=""
                    ),
                    collapse=""
                )
                
                # summary by element
                return(paste(paste(x,elems,"\n",sep =""),collapse="/n"))
            }
        },"")

        # cat some text
        cat("- object class: enrich_GO_terms",
            "\n- ontology: ",slot(object,"ont"),
            "\n- method: ",slot(object,"method"),
            "\n- summary:\n ", obj_summary,
            "- enrichment pvalue cutoff:",
            paste("\n       ",Data$conditions,":",slot(object,"cutoff")[[1]]),
            "\n- enrich GOs (in at least one list): ",nrow(slot(object,"data"))," GO terms of ",nrow(Data)," conditions.",
            paste("\n       ",Data$conditions,":",Data$`significant GO terms number`,"terms"),sep=""
        )
    }
)
