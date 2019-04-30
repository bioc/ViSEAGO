#' @title GO_clusters class object
#' @description This class is invoked by \code{\link{GOterms_heatmap}} and  \code{\link{GOclusters_heatmap}} methods to store all results produced.
#' @family GO_clusters
#' @slot db database source.
#' @slot stamp date of stamp.
#' @slot organism target species.
#' @slot topGO topGO objects summary.
#' @slot ont ontology (MF, BP, or CC).
#' @slot IC Information Content (IC).
#' @slot enrich_GOs \code{\link{enrich_GO_terms-class}} object.
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
        db="character",
        stamp = "character",
        organism="character",
        ont="character",
        topGO="list",
        IC="numeric",
        enrich_GOs="enrich_GO_terms",
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

        # Extract table
        Data<-slot(
            slot(
                object,
                "enrich_GOs"
            ),
        "data"
        )

        # Extract pvalues
        Data<-Data[,grep("\\.pvalue",names(Data)),with=FALSE]

        # count significant pvalues by condition
        Data<-Data[,lapply(.SD,function(x){sum(x<0.01,na.rm = T)}),.SDcols=seq_len(ncol(Data))]

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
        cat("- object class: GO_clusters",
            "\n- database: ",slot(object,"db"),
            "\n- stamp/version: ",slot(object,"stamp"),
            "\n- organism id: ",slot(object,"organism"),
            "\n- ontology: ",slot(object,"ont"),
            "\n- input:\n        ",
            paste(
                paste(
                    names(slot(slot(object,"enrich_GOs"),"input")),
                    vapply(slot(slot(object,"enrich_GOs"),"input"),function(x){paste(x,collapse=", ")},""),
                    sep=": "
                ),
                collapse="\n        "
            ),
            "\n- topGO summary:\n ",
            topGO,
            "\n- enrich GOs data.table: ",
            nrow(slot(slot(object,"enrich_GOs"),"data")),
            " GO terms of ",
            length(
                grep(
                    "\\.pvalue",
                    names(slot(slot(object,"enrich_GOs"),"data"))
                )
            ),
            " conditions.",
            paste(
                "\n       ",
                Data$conditions,
                ":",
                 Data$`significant GO terms number`,
                "terms"
                ),
                "\n- clusters distances: ",
                paste(
                    names(
                        slot(object,"clusters_dist")
                    ),
                    collapse=", "
                ),
                "\n- Heatmap:",
                "\n          * GOterms: ",
                !is.null(
                    slot(object,"heatmap")$GOterms
                ),
                "\n                    - GO.tree:\n                              ",
                paste(
                    paste(
                        names(
                            unlist(
                                slot(object,"hcl_params")$GO.tree
                            )
                        ),
                        unlist(
                            slot(object,"hcl_params")$GO.tree
                        ),sep=": "
                    ),
                    collapse="\n                              "
                ),
                "\n                              number of clusters: ",
                length(
                    unique(
                        slot(slot(object,"enrich_GOs"),"data")[,"GO.cluster",with=FALSE]
                    )
                ),
                "\n                              clusters min size: ",
                round(
                    min(
                        slot(slot(object,"enrich_GOs"),"data")[,"GO.cluster",with=FALSE],
                        na.rm=TRUE
                    ),
                    digits=0
                ),
                "\n                              clusters mean size: ",
                round(
                    mean(
                        unlist(slot(slot(object,"enrich_GOs"),"data")[,"GO.cluster",with=FALSE]),
                        na.rm=TRUE
                    ),
                    digits=0
                ),
                "\n                              clusters max size: ",
                round(
                    max(
                        slot(slot(object,"enrich_GOs"),"data")[,"GO.cluster",with=FALSE],
                        na.rm=TRUE
                    ),
                    digits=0
                ),
                "\n                   - sample.tree: ",
                paste(
                    paste(
                        names(
                            unlist(
                                slot(object,"hcl_params")$samples.tree
                            )
                        ),
                        unlist(
                            slot(object,"hcl_params")$samples.tree
                        ),
                        sep=": "
                    ),
                    collapse="\n                                 "
                ),
                if(is.null(slot(object,"hcl_params")$samples.tree)){"FALSE"},
                "\n          * GOclusters: ",
                !is.null(
                    slot(object,"heatmap")$GOclusters
                ),
                if(!is.null(slot(object,"heatmap")$GOclusters)){
                    paste(
                        "\n                       - tree:\n                             ",
                        paste(
                            paste(
                                names(
                                    unlist(
                                        slot(object,"hcl_params")$GO.clusters
                                    )
                                ),
                                unlist(
                                    slot(object,"hcl_params")$GO.clusters
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
