#' @title Display an interactive or static heatmap.
#' @description Display a heatmap in interactive or static mode.
#' @importFrom methods setGeneric setMethod slot is
#' @importFrom plotly export
#' @importFrom webshot install_phantomjs
#' @family enrich_GO_terms
#' @family GO_clusters
#' @family visualization
#' @param object a \code{\link{GO_clusters-class}} object from \code{\link{GOterms_heatmap}} or \code{\link{GOclusters_heatmap}}.
#' @param type could be "GOterms" to display GOterms clustering heatmap, or "GOclusters" to display GOclusters heatmap.
#' @param file static image output file name (default to NULL).
#' @details
#' This method displays an interactive heatmap (if \code{file}=NULL) from \code{\link{GO_clusters-class}} object for "GOterms" or "GOclusters" type.\cr
#' A static png image could be printed by setting \code{file} argument.
#' @examples
#' ##################
#' # Display GO terms heatmap
#' ViSEAGO::show_heatmap(Wang_clusters_wardD2,"GOterms")
#'
#' ##################
#' # Display GO clusters heatmap
#' ViSEAGO::show_heatmap(Wang_clusters_wardD2,"GOclusters")
#'
#' ##################
#' # Print GO terms heatmap
#' ViSEAGO::show_heatmap(Wang_clusters_wardD2,"GOterms","./data/output/GOterms_heatmap.png")
#'
#' ##################
#' # Print GO clusters heatmap
#' ViSEAGO::show_heatmap(Wang_clusters_wardD2,"GOclusters","./data/output/GOclusters_heatmap.png")
#' @export
setGeneric(name="show_heatmap",def=function(object,type,file=NULL) {standardGeneric("show_heatmap")})
#' @importFrom methods setMethod
setMethod("show_heatmap",signature="GO_clusters",definition=function(object,type,file=NULL) {

  ##################
  # check type argument
  type=base::match.arg(type,c("GOterms","GOclusters"))

  ##################
  # heatmap according type
  heatmap<-base::switch(type,

    ##################
    # GOTermsHeatmap
    GOterms=methods::slot(object,"heatmap")$GOterms,

    ##################
    # GOclusters_heatmap
    GOclusters=methods::slot(object,"heatmap")$GOclusters
  )

  ##################
  # return interactive or static heatmap according file
  if(base::is.null(file)){

    ##################
    # display heatmap
    heatmap

  }else{

    if(type=="GOterms"){

      ##################
      # number of rows
      rowlen=base::nrow(
        methods::slot(
          methods::slot(
            object,
            "enrich_GOs"
          ),
          "data"
        )
      )

      ##################
      # adjust minimum size
      if(rowlen<10){rowlen=10}

      ##################
      # compute height
      rowlen=rowlen^(1.70+1.70*exp(-rowlen/20))

      ##################
      # max height limit
      if(rowlen>10000){rowlen=10000}

      ##################
      # adjust heatmpa size
      heatmap<-plotly::layout(
        heatmap,
        height=rowlen
      )
    }

    ##################
    # print heatmap
    plotly::export(heatmap,file=file)
  }
})
