#' @title GO_clusters class object
#' @description This class is invoked by \code{\link{GOterms_heatmap}} and  \code{\link{GOclusters_heatmap}} methods to store all results produced.
#' @importFrom methods setClass slot
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
setClass("GO_clusters",
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
#' @importFrom methods setMethod
setMethod("show", "GO_clusters",function(object) {

  ###################
  # Extract table
  Data<-methods::slot(
    methods::slot(
      object,
      "enrich_GOs"
    ),
  "data"
  )

  ###################
  # Extract pvalues
  Data<-Data[,base::grep("\\.pvalue",base::names(Data)),with=F]

  ###################
  # count significant pvalues by condition
  Data<-Data[,lapply(.SD,function(x){base::sum(x<0.01,na.rm = T)}),
             .SDcols=1:base::ncol(Data)
             ]

  ###################
  # melt the table
  Data<-data.table::melt.data.table(
    Data,
    measure.vars=base::names(Data),
    variable.name = "conditions",
    value.name = "significant GO terms number"
  )

  ###################
  # remove .pvalue in conditions column
  Data[,conditions:=base::gsub("\\.pvalue","",conditions)]


  ###################
  # get topGO information
  topGO<-methods::slot(object,"topGO")

  ###################
  # format topGO information
  topGO<-base::sapply(base::names(topGO),function(x){

    ###################
    # topGO subset
    Data<-topGO[[x]]

    ###################
    # summery by element
    elems<-base::paste(base::sapply(base::names(Data),function(y){

      ###################
      # summery by element
      base::paste(
        base::paste("   ",y),
        "\n       ",
        base::paste(
          base::paste(
            base::names(Data[[y]]),
            base::sub("\n.+$","",base::unlist(Data[[y]])),
            sep=": "),
          collapse="\n        "
        )
      )
    }),collapse="\n  ")

    ###################
    # summery by element
    base::paste(base::paste(x,elems,sep ="\n  "),"\n ")
  })

  ###################
  # cat some text
  base::cat("- object class: GO_clusters",
    "\n- database: ",methods::slot(object,"db"),
    "\n- stamp/version: ",methods::slot(object,"stamp"),
    "\n- organism id: ",methods::slot(object,"organism"),
    "\n- ontology: ",methods::slot(object,"ont"),
    "\n- input:\n        ",
      paste(
        paste(base::names(
          methods::slot(
            methods::slot(
              object,
              "enrich_GOs"
            ),
          "input")
        ),
        base::sapply(
          methods::slot(
            methods::slot(
              object,
              "enrich_GOs"
              ),
            "input"
          ),
          function(x){base::paste(x,collapse=", ")}
        ),
        sep=": "),
        collapse="\n        "
      ),
    "\n- topGO summary:\n ",
    topGO,
    "\n- enrich GOs data.table: ",
      base::nrow(
        methods::slot(
          methods::slot(
            object,
            "enrich_GOs"
          ),
          "data"
        )
      ),
    " GO terms of ",
    base::length(
      base::grep(
        "\\.pvalue",
        base::names(
          methods::slot(
            methods::slot(
              object,
              "enrich_GOs"
            ),
            "data"
          )
        )
      )
    ),
    " conditions.",
    base::paste(
      "\n       ",
      Data$conditions,
      ":",
      Data$`significant GO terms number`,
      "terms"
    ),
    "\n- clusters distances: ",
    base::paste(
      base::names(
        methods::slot(object,"clusters_dist")
      ),
      collapse=", "
    ),
    "\n- Heatmap:",
    "\n          * GOterms: ",
    !base::is.null(
      methods::slot(object,"heatmap")$GOterms
    ),
    "\n                    - GO.tree:\n                              ",
    paste(
      base::paste(
        base::names(
          base::unlist(
            methods::slot(object,"hcl_params")$GO.tree
          )
        ),
        base::unlist(
          methods::slot(object,"hcl_params")$GO.tree
        ),sep=": "),
      collapse="\n                              "
    ),
    "\n                              number of clusters: ",
    base::length(
      base::unique(
        methods::slot(
          methods::slot(
            object,
            "enrich_GOs"
          ),
          "data"
        )[,GO.cluster]
      )
    ),
    "\n                              clusters min size: ",
    base::round(
      base::min(
        methods::slot(
          methods::slot(
            object,
            "enrich_GOs"
          ),
          "data"
        )[,GO.cluster]
      ),
      digits=0
    ),
    "\n                              clusters mean size: ",
    base::round(
      base::mean(
        methods::slot(
          methods::slot(
            object,
            "enrich_GOs"
          ),
          "data"
        )[,GO.cluster]
      ),
      digits=0
    ),
    "\n                              clusters max size: ",
    base::round(
      base::max(
        methods::slot(
          methods::slot(
            object,
            "enrich_GOs"
          ),
          "data"
        )[,GO.cluster]
      ),
      digits=0
    ),
    "\n                   - sample.tree: ",
    paste(
      base::paste(
        base::names(
          base::unlist(
            methods::slot(object,"hcl_params")$samples.tree
          )
        ),
        base::unlist(
          object@hcl_params$samples.tree
        ),
        sep=": "
      ),collapse="\n                                 "
    ),
    if(base::is.null(methods::slot(object,"hcl_params")$samples.tree)){"FALSE"},
    "\n          * GOclusters: ",
    !base::is.null(
      methods::slot(object,"heatmap")$GOclusters
    ),
    if(!base::is.null(methods::slot(object,"heatmap")$GOclusters)){
      paste(
        "\n                       - tree:\n                             ",
        base::paste(
          base::paste(
            base::names(
              base::unlist(
                methods::slot(object,"hcl_params")$GO.clusters
              )
            ),
          base::unlist(
            methods::slot(object,"hcl_params")$GO.clusters
          ),
          sep=": "
        ),
        collapse="\n                              "),
      collapse=""
    )
    },
    sep=""
  )
})
