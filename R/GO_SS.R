#' @title GO_SS class object definition.
#' @description This class is invoked by \code{\link{build_GO_SS}} method in order to store \code{\link{enrich_GO_terms-class}} object, Information Content (IC),
#' and GO terms or groups distances objects based on semantic similarity.
#' @family GO_semantic_similarity
#' @slot ont ontology used "MF", "BP", or "CC".
#' @slot enrich_GOs \code{\link{merge_enrich_terms}} output object (\code{\link{enrich_GO_terms-class}} object).
#' @slot IC Information Content (IC)
#' @slot terms_dist \code{list} of GO terms or groups distances objects based on semantic similarity.
#' @include enrich_GO_terms.R
setClass(
    "GO_SS",
    slots=c(
        ont="character",
        enrich_GOs="enrich_GO_terms",
        IC ="numeric",
        terms_dist="list"
    )
)

#' @aliases GO_SS
#' @importFrom data.table melt.data.table .SD
setMethod(
    "show",
    signature="GO_SS",
    function(object){

        # keep full object
        full_object<-object
        
        # extract object
        object<-slot(object,"enrich_GOs")

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
                return(paste(paste("\n ",x,elems,sep =""),collapse="/n"))
            }
        },"")

        # cat some text
        cat("- object class: GO_SS",
            "\n- ontology: ",slot(object,"ont"),
            "\n- method: ",slot(object,"method"),
            "\n- summary:\n", obj_summary,
            "- enrichment pvalue cutoff:",
            paste("\n       ",Data$conditions,":",slot(object,"cutoff")[[1]]),
            "\n- enrich GOs (in at least one list): ",nrow(slot(object,"data"))," GO terms of ",nrow(Data)," conditions.",
            paste("\n       ",Data$conditions,":",Data$`significant GO terms number`,"terms"),
            if(length(slot(full_object,"terms_dist"))>0){
                paste(
                    "\n- terms distances: ",
                    paste(
                        names(
                            slot(full_object,"terms_dist")
                        ),
                        collapse=", "
                    )
                )
            },
            sep=""
        )
    }
)
