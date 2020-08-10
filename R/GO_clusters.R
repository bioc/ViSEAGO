#' @title GO_clusters class object
#' @description This class is invoked by \code{\link{GOterms_heatmap}} and  \code{\link{GOclusters_heatmap}} methods to store all results produced.
#' @family GO_clusters
#' @slot ont ontology used "MF", "BP", or "CC".
#' @slot enrich_GOs \code{\link{enrich_GO_terms-class}} object.
#' @slot IC Information Content (IC).
#' @slot terms_dist distance between GO terms based on semantic similiarity.
#' @slot clusters_dist distance between GO groups based on semantic similiarity.
#' @slot hcl_params Hierarchical clustering parameters used.
#' @slot dendrograms GO terms and samples \code{dendrograms}.
#' @slot samples.gp samples groups.
#' @slot heatmap GO terms and GO groups heatmaps.
#' @include enrich_GO_terms.R
setClass(
    "GO_clusters",
    slots=c(
        ont="character",
        enrich_GOs="enrich_GO_terms",
        IC ="numeric",
        terms_dist="list",
        clusters_dist="list",
        hcl_params="list",
        dendrograms="list",
        samples.gp="numeric",
        heatmap="list"
    )
)

#' @aliases GO_clusters
#' @importFrom data.table .SD melt.data.table :=
setMethod(
    "show",
    signature="GO_clusters",
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
        cat("- object class: GO_clusters",
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
            "\n- clusters distances: ",
            paste(
                names(
                    slot(full_object,"clusters_dist")
                ),
                collapse=", "
            ),
            "\n- Heatmap:",
            "\n          * GOterms: ",
            !is.null(
                slot(full_object,"heatmap")$GOterms
            ),
            "\n                    - GO.tree:\n                              ",
            paste(
                paste(
                    names(
                        unlist(
                            slot(full_object,"hcl_params")$GO.tree
                        )
                    ),
                    unlist(
                        slot(full_object,"hcl_params")$GO.tree
                    ),sep=": "
                ),
                collapse="\n                              "
            ),
            "\n                              number of clusters: ",
            length(
                unique(
                    unlist(slot(slot(full_object,"enrich_GOs"),"data")[,"GO.cluster",with=FALSE])
                )
            ),
            "\n                              clusters min size: ",
            round(
                min(
                    table(slot(slot(full_object,"enrich_GOs"),"data")[,"GO.cluster",with=FALSE]),
                    na.rm=TRUE
                ),
                digits=0
            ),
            "\n                              clusters mean size: ",
            round(
                mean(
                    table(slot(slot(full_object,"enrich_GOs"),"data")[,"GO.cluster",with=FALSE]),
                    na.rm=TRUE
                ),
                digits=0
            ),
            "\n                              clusters max size: ",
            round(
                max(
                    table(slot(slot(full_object,"enrich_GOs"),"data")[,"GO.cluster",with=FALSE]),
                    na.rm=TRUE
                ),
                digits=0
            ),
            "\n                   - sample.tree: ",
            paste(
                paste(
                    names(
                        unlist(
                            slot(full_object,"hcl_params")$samples.tree
                        )
                    ),
                    unlist(
                        slot(full_object,"hcl_params")$samples.tree
                    ),
                    sep=": "
                ),
                collapse="\n                                 "
            ),
            if(is.null(slot(full_object,"hcl_params")$samples.tree)){"FALSE"},
            "\n          * GOclusters: ",
            !is.null(
                slot(full_object,"heatmap")$GOclusters
            ),
            if(!is.null(slot(full_object,"heatmap")$GOclusters)){
                paste(
                    "\n                       - tree:\n                             ",
                    paste(
                        paste(
                            names(
                                unlist(
                                    slot(full_object,"hcl_params")$GO.clusters
                                )
                            ),
                            unlist(
                                slot(full_object,"hcl_params")$GO.clusters
                            ),
                            sep=": "
                        ),
                        collapse="\n                              "
                    ),
                    collapse=""
                )
            },
            sep=""
        )
    }
)
