#' @title Build a clustering heatmap on GO groups.
#' @description This method computes a clustering heatmap based on GO groups semantic similarity.
#' @importFrom data.table data.table .N
#' @importFrom methods setGeneric setMethod signature
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom plotly layout
#' @importFrom heatmaply heatmaply
#' @importFrom dendextend branches_attr_by_clusters set get_leaves_attr
#' @importFrom RColorBrewer brewer.pal
#' @family GO_clusters
#' @family semantic_similiarity
#' @family visualization
#' @param object a  \code{\link{GO_clusters-class}} object from \code{\link{compute_SS_distances}}.
#' @param tree a named list with:
#'  \describe{
#'      \item{distance ("BMA" by default)}{
#'      distance computed from the semantic similarity for GO groups which could be "max", "avg", "rcmax",or "BMA".}
#'      \item{aggreg.method ("ward.D2" by default)}{aggregation method criteria from \code{\link[stats]{hclust}} (ward.D",
#'      "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid") to build a \code{dendrogram}.}
#'      \item{rotate}{sort the branches of the tree based on a vector - eithor of labels order or the labels in their new order}
#'  }
#' @details This method computes a clustering heatmap based on GO groups semantic similarity (computed with \code{\link{compute_SS_distances}}).\cr
#' The heatmap color intensity corresponds to the number of GO terms in each GO group.\cr
#' GO group description is defined as the first common GO ancestor with the cluster identifiant in brackets.\cr
#' The dendrogram branches are colored according to GO terms clusters.
#' @references
#' Matt Dowle and Arun Srinivasan (2017). data.table: Extension of `data.frame`. R package version 1.10.4. https://CRAN.R-project.org/package=data.table.
#'
#' Tal Galili (2015). dendextend: an R package for visualizing, adjusting, and comparing trees of hierarchical clustering.
#' Bioinformatics. DOI:10.1093/bioinformatics/btv428.
#'
#' Tal Galili (2017). heatmaply: Interactive Cluster Heat Maps Using 'plotly'.
#' R package version 0.9.1. https://CRAN.R-project.org/package=heatmaply.
#'
#' Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1-2. https://CRAN.R-project.org/package=RColorBrewer.
#'
#' Carson Sievert, Chris Parmer, Toby Hocking, Scott Chamberlain, Karthik Ram, Marianne Corvellec and Pedro Despouy (2017).
#' plotly: Create Interactive Web Graphics via 'plotly.js'. R package version 4.6.0. https://CRAN.R-project.org/package=plotly.
#'
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
#' @return a \code{\link{GO_clusters-class}} object.
#' @include GO_clusters.R
#' @examples
#' ###################
#' # load data example
#' utils::data(
#'  myGOs,
#'  package="ViSEAGO"
#' )
#' \dontrun{
#' ###################
#' # compute GO terms Semantic Similarity distances
#' myGOs<-ViSEAGO::compute_SS_distances(
#'  myGOs,
#'  distance="Wang"
#' )
#'
#' ##################
#' # GOtermsHeatmap with default parameters
#' Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#'  myGOs,
#'  showIC=TRUE,
#'  showGOlabels=TRUE,
#'  GO.tree=base::list(
#'   tree=base::list(
#'    distance="Wang",
#'    aggreg.method="ward.D2",
#'    rotate=NULL
#'   ),
#'   cut=base::list(
#'    dynamic=base::list(
#'     pamStage=TRUE,
#'     pamRespectsDendro=TRUE,
#'     deepSplit=2,
#'     minClusterSize =2
#'    )
#'   )
#'  ),
#'  samples.tree=NULL
#' )
#'
#' ###################
#' # compute clusters of GO terms Semantic Similarity distances
#' Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(
#'  Wang_clusters_wardD2,
#'  distance="BMA"
#' )
#'
#' ###################
#' # GOclusters heatmap
#' Wang_clusters_wardD2<-ViSEAGO::GOclusters_heatmap(
#'  Wang_clusters_wardD2,
#'  tree=list(
#'   distance="BMA",
#'   aggreg.method="ward.D2",
#'   rotate=NULL
#'  )
#' )
#' }
#' @exportMethod GOclusters_heatmap
#' @name GOclusters_heatmap
#' @rdname GOclusters_heatmap-methods
#' @exportMethod GOclusters_heatmap
setGeneric(
  name="GOclusters_heatmap",
  def=function(object,tree=base::list(distance="BMA",aggreg.method="ward.D2",rotate=NULL)){
    base::standardGeneric("GOclusters_heatmap")
  }
)

#' @rdname GOclusters_heatmap-methods
#' @aliases GOclusters_heatmap
setMethod("GOclusters_heatmap",
  methods::signature(
    object="GO_clusters",
    tree="list"
  ),
  definition=function(object,tree){

  #################
  # if empty
  if(is.null(methods::slot(object,"clusters_dist"))){
    stop("Please use GO_clusters class object with computed clusters distances using ViSEAGO::compute_SS_distance()")
  }

  #################
  # if empty
  if(!tree$distance%in%base::names(methods::slot(object,"clusters_dist"))){
    stop(paste("Please use GO_clusters class object with",tree$distance,"computed clusters distances using ViSEAGO::compute_SS_distance()"))
  }

  ###################
  # add calculated distance
  dist<-methods::slot(object,"clusters_dist")[[tree$distance]]

  ###################
  # perform hclust
  Tree<-stats::hclust(
    dist,
    method =tree$aggreg.method
  )

  ###################
  # create dendrogram
  dend<-stats::as.dendrogram(Tree)

  ###################
  # rotate
  if(!base::is.null(tree$rotate)){
    dend<-dendextend::rotate(
      dend,
      tree$rotate
    )
  }

  ###################
  # create dendrogram
  ord<-stats::order.dendrogram(dend)

  ###################
  # color labels of cutting tree
  colors= base::unique(
    dendextend::get_leaves_attr(
      methods::slot(object,"dendrograms")$GO,
      "edgePar"
    )
  )

  ###################
  # extract GO terms and clusters
  mat<-methods::slot(
    methods::slot(
      object,
      "enrich_GOs"
      ),
    "data"
  )[,base::c("GO.cluster","GO.ID"),with=FALSE]

  ###################
  # count terms by clusters
  mat<-mat[,base::list("count"=.N),by="GO.cluster"]

  ###################
  # count terms by clusters
  mat<-base::as.matrix(mat[,"count",with=FALSE])

  ###################
  # colors table
  colors=data.table::data.table(
    gp=base::seq_len(base::nrow(mat)),
    color=colors
  )

  ###################
  # merge cluster term assignation and corresponding color
  colors<-merge(
    data.table::data.table(gp=colors$gp[ord]),
    colors,
    by="gp",
    all.x=TRUE,
    sort=FALSE
  )

  ###################
  # color branches according clusters
  dend<-dendextend::branches_attr_by_clusters(
    dend,
    base::seq_len(base::nrow(mat)),
    values =colors$color
  )

  ###################
  # assign text color
  dend<-dendextend::set(dend,"labels_col",colors$color)

  ###################
  # create dendrogram
  dend<-dendextend::set(dend,"labels_cex",0.3)

  ###################
  # extract GO terms and clusters
  base::row.names(mat)<-base::attr(dist,"Labels")

  ###################
  # custom row names
  base::row.names(mat)<-base::paste("<br>cluster:", base::row.names(mat))
  base::row.names(mat)<-base::gsub("_GO:","<br>GO.ID: GO:",base::row.names(mat))
  base::row.names(mat)<-base::gsub("_","<br>GO.name: ", base::row.names(mat))

  ###################
  # draw heatmapply
  hm<-heatmaply::heatmaply(

    ###################
    # the initial matrix
    x=mat,

    ###################
    # row label
    labRow=base::row.names(mat),

    ###################
    # row label
    labCol=base::colnames(mat),

    ###################
    # the row dendrogram
    Rowv=dend,

    ###################
    # the ordered matrix according dendrograms for columns
    Colv=FALSE,

    ###################
    # the color palette
    scale_fill_gradient_fun=ggplot2::scale_fill_gradient(
      name="GO terms count",
      low="#FCFBFD",
      high="#3F007D"
    ),

    ###################
    # the width of dendrogramm
    branches_lwd = 0.4,

    ###################
    # color bar length
    colorbar_len=0.05
  )

  ###################
  # remove row tag
  hm$x$data[[1]]$text<-base::gsub("<br>row: ","",hm$x$data[[1]]$text)

  ###################
  # custom row text
  row.text=base::gsub("^.+GO.name: ","",base::rev(base::row.names(mat)[ord]))

  ###################
  # cut very long definition
  row.text[nchar(row.text)>50]<-base::paste(base::substring(
  row.text[nchar(row.text)>50],1,50),"...",sep="")

  ###################
  # add cluster number in brackets
  row.text=base::paste(row.text," (cl",
  base::gsub("^<br>cluster: |<br>GO.ID.+$","",
  base::rev(base::row.names(mat)[ord])),")",sep="")

  hm<-plotly::layout(hm,

    #################
    # add title
    title=paste(tree$distance,"GOclusters distance heatmap"),

    #################
    # x axis
    xaxis=base::list(
      domain=c(0,0.2),
      family="Times New Roman",
      tickfont=base::list(size=10)
    ),
    xaxis2=base::list(
      domain=c(0.2,0.8)
    ),

    #################
    # title size
    font=base::list(size=14),

    #################
    # y axis
    yaxis=base::list(
      family="Times New Roman",
      tickmode="array",
      tickvals=base::seq_len(base::nrow(mat)),
      ticktext=row.text,
      tickfont=base::list(size=10)
    ),
    margin = list(l =300,r=0, b =150,t=50)
  )

  ###################
  #  hm to list
  hm<-base::list(hm)

  ###################
  # give names to hm list
  base::names(hm)<-"GOclusters"

  ###################
  # give names to hm list
  methods::slot(object,"heatmap")<-c(methods::slot(object,"heatmap"),hm)

  ###################
  # add hcl params
  methods::slot(object,"hcl_params")<-base::c(methods::slot(object,"hcl_params"),GO.clusters=base::list(tree))

  ###################
  # return the object
  base::return(object)
})
