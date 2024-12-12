#' @title Display an interactive or static heatmap.
#' @description Display a heatmap in interactive or static mode.
#' @importFrom plotly export layout
#' @importFrom grDevices dev.off png
#' @importFrom graphics text
#' @importFrom stats end start
#' @family enrich_GO_terms
#' @family GO_clusters
#' @family visualization
#' @param object a \code{\link{GO_clusters-class}} object from \code{\link{GOterms_heatmap}} or \code{\link{GOclusters_heatmap}}.
#' @param type could be "GOterms" to display GOterms clustering heatmap, or "GOclusters" to display GOclusters heatmap.
#' @param file static png output file name (default to NULL).
#' @param width static image width (default to 1000).
#' @param height static image height (default to 1000).
#' @param plotly_update update plotly html dependencies (default to FALSE).
#' @details
#' This method displays an interactive heatmap (if \code{file}=NULL) from \code{\link{GO_clusters-class}} object for "GOterms" or "GOclusters" type.\cr
#' A static png image could be printed by setting \code{file} argument.\cr
#' Interactive heatmap cannot be displayed between two R versions.
#' Then interactive view (build with previous R version) can be updated to new R version using \code{plotly_update} argument setting to TRUE.
#' @return display or print heatmap.
#' @examples
#' # load data example
#' data(
#'     myGOs,
#'     package="ViSEAGO"
#' )
#' \dontrun{
#' # compute GO terms Semantic Similarity distances
#' myGOs<-ViSEAGO::compute_SS_distances(
#'     myGOs,
#'     distance="Wang"
#' )
#'
#' # build MDS plot for a GO_SS-class distance object
#' ViSEAGO::MDSplot(myGOs)
#'
#' # GOtermsHeatmap with default parameters
#' Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#'     myGOs,
#'     showIC=TRUE,
#'     showGOlabels=TRUE,
#'     GO.tree=list(
#'         tree=list(
#'             distance="Wang",
#'             aggreg.method="ward.D2",
#'             rotate=NULL
#'         ),
#'         cut=list(
#'             dynamic=list(
#'                 pamStage=TRUE,
#'                 pamRespectsDendro=TRUE,
#'                 deepSplit=2,
#'                 minClusterSize =2
#'             )
#'         )
#'     ),
#'     samples.tree=NULL
#' )
#'
#' # Display GO terms heatmap
#' ViSEAGO::show_heatmap(
#'     Wang_clusters_wardD2,
#'     "GOterms"
#' )
#'
#' # Print GO terms heatmap
#' ViSEAGO::show_heatmap(
#'     Wang_clusters_wardD2,
#'     "GOterms",
#'     "GOterms_heatmap.png"
#' )
#'
#' # compute clusters of GO terms Semantic Similarity distances
#' Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(
#'     Wang_clusters_wardD2,
#'     distance="BMA"
#' )
#'
#' # GOclusters heatmap
#' Wang_clusters_wardD2<-ViSEAGO::GOclusters_heatmap(
#'     Wang_clusters_wardD2,
#'     tree=list(
#'         distance="BMA",
#'         aggreg.method="ward.D2",
#'         rotate=NULL
#'     )
#' )
#'
#' # Display GO clusters heatmap
#' ViSEAGO::show_heatmap(
#'     Wang_clusters_wardD2,
#'     "GOclusters"
#' )
#'
#' # Print GO clusters heatmap
#' ViSEAGO::show_heatmap(
#'     Wang_clusters_wardD2,
#'     "GOclusters",
#'     "GOclusters_heatmap.png"
#' )
#' }
#' @name show_heatmap
#' @rdname show_heatmap-methods
#' @exportMethod show_heatmap
setGeneric(
    name="show_heatmap",
    def=function(object,type,file=NULL,plotly_update=FALSE,height=1000,width=600) {
        standardGeneric("show_heatmap")
    }
)

#' @rdname show_heatmap-methods
#' @aliases show_heatmap
setMethod(
    "show_heatmap",
    signature(
        object="GO_clusters",
        type="character"
    ),
    definition=function(object,type,file,plotly_update,height,width){

        # check type argument
        type=match.arg(type,c("GOterms","GOclusters"))

        # heatmap according type
        heatmap<-switch(type,

            # GOTermsHeatmap
            GOterms=slot(object,"heatmap")$GOterms,

            # GOclusters_heatmap
            GOclusters=slot(object,"heatmap")$GOclusters
        )

        # redirect plotly html_dependency to local packages
        if(plotly_update==TRUE){

            for(i in seq_len(length(heatmap[["dependencies"]]))){
                heatmap[["dependencies"]][[i]][["src"]][["file"]]<-gsub(
                    "^.+[[:digit:]]\\.[[:digit:]]",
                    .libPaths()[1],
                    heatmap[["dependencies"]][[i]][["src"]][["file"]]
                )
            }
        }

        # return interactive or static heatmap according file
        if(is.null(file)){

            # display heatmap
            heatmap

        }else{

            # static heatmap according type
            heatmap<-switch(
                type,
                GOterms=slot(object,"heatmap")$GOterms_static,
                GOclusters=slot(object,"heatmap")$GOclusters_static
            )

            # plot
            png(file,width=width,height=height)
                print(
                    heatmap
                )
            dev.off()
        }
    }
)
