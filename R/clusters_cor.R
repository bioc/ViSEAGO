#' @title Compute distance matrix between dendrograms partitions.
#' @description Build a distance or correlation matrix between partitions from dendrograms.
#' @importFrom data.table data.table
#' @importFrom igraph compare
#' @importFrom methods setGeneric setMethod signature
#' @family GO_clusters
#' @param clusters  a \code{list} of \code{\link{GO_clusters-class}} objects,
#' from \code{\link{GOterms_heatmap}} or \code{\link{GOclusters_heatmap}}, named as \code{character}.
#' @param method available methods  ("vi", "nmi", "split.join", "rand", or "adjusted.rand") from \pkg{igraph} package \code{\link[igraph]{compare}} function.
#' @return a distance or a correlation matrix.
#' @references
#' Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. http://igraph.org.
#' @include GO_clusters.R
#' @examples
#' ###################
#' # load example object
#' utils::data(
#'  myGOs,
#'  package="ViSEAGO"
#' )
#'
#' \dontrun{
#' ###################
#' # compute Semantic Similarity (SS)
#' myGOs<-ViSEAGO::compute_SS_distances(
#'  myGOs,
#'  distance=base::c("Resnik","Lin","Rel","Jiang","Wang")
#' )
#'
#' ##################
#' # Resnik distance GO terms heatmap
#' Resnik_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#'   myGOs,
#'   showIC=TRUE,
#'   showGOlabels=TRUE,
#'   GO.tree=base::list(
#'     tree=base::list(
#'       distance="Resnik",
#'       aggreg.method="ward.D2"
#'     ),
#'     cut=base::list(
#'       dynamic=base::list(
#'         deepSplit=2,
#'         minClusterSize =2
#'       )
#'     )
#'   ),
#'   samples.tree=NULL
#' )
#'
#' ##################
#' # Lin distance GO terms heatmap
#' Lin_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#'   myGOs,
#'   showIC=TRUE,
#'   showGOlabels=TRUE,
#'   GO.tree=base::list(
#'     tree=base::list(
#'       distance="Lin",
#'       aggreg.method="ward.D2"
#'     ),
#'     cut=base::list(
#'       dynamic=base::list(
#'         deepSplit=2,
#'         minClusterSize =2
#'       )
#'     )
#'   ),
#'   samples.tree=NULL
#' )
#'
#' ##################
#' # Rel distance GO terms heatmap
#' Rel_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#'   myGOs,
#'   showIC=TRUE,
#'   showGOlabels=TRUE,
#'   GO.tree=base::list(
#'     tree=base::list(
#'       distance="Rel",
#'       aggreg.method="ward.D2"
#'     ),
#'     cut=base::list(
#'       dynamic=base::list(
#'         deepSplit=2,
#'         minClusterSize =2
#'       )
#'     )
#'   ),
#'   samples.tree=NULL
#' )
#'
#' ##################
#' # Jiang distance GO terms heatmap
#' Jiang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#'   myGOs,
#'   showIC=TRUE,
#'   showGOlabels=TRUE,
#'   GO.tree=base::list(
#'     tree=base::list(
#'       distance="Jiang",
#'       aggreg.method="ward.D2"
#'     ),
#'     cut=base::list(
#'       dynamic=base::list(
#'         deepSplit=2,
#'         minClusterSize =2
#'       )
#'     )
#'   ),
#'   samples.tree=NULL
#' )
#'
#' ##################
#' # Wang distance GO terms heatmap
#' Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#'   myGOs,
#'   showIC=TRUE,
#'   showGOlabels=TRUE,
#'   GO.tree=base::list(
#'     tree=base::list(
#'       distance="Wang",
#'       aggreg.method="ward.D2"
#'     ),
#'     cut=base::list(
#'       dynamic=base::list(
#'         deepSplit=2,
#'         minClusterSize =2
#'       )
#'     )
#'   ),
#'   samples.tree=NULL
#' )
#' }
#'
#' ###################
#' # clusters to compare
#' clusters<-base::list(
#'  Resnik="Resnik_clusters_wardD2",
#'  Lin="Lin_clusters_wardD2",
#'  Rel="Rel_clusters_wardD2",
#'  Jiang="Jiang_clusters_wardD2",
#'  Wang="Wang_clusters_wardD2"
#' )
#'
#' \dontrun{
#' ###################
#' # global dendrogram clustering correlation
#' clust_cor<-ViSEAGO::clusters_cor(
#'  clusters,
#'  method="adjusted.rand"
#' )
#' }
#' @name clusters_cor
#' @rdname clusters_cor-methods
#' @exportMethod clusters_cor
setGeneric(name="clusters_cor",def=function(clusters,method="adjusted.rand") {standardGeneric("clusters_cor")})

#' @rdname clusters_cor-methods
#' @aliases clusters_cor
setMethod("clusters_cor",
  methods::signature(
    clusters="list",
    method="character"
  ),
  definition=function(clusters,method) {

  ###################
  # method
  method=base::match.arg(method,c("vi", "nmi", "split.join", "rand","adjusted.rand"))

  ###################
  # extract clusters and GO.ID
  ###################

    ###################
    # extract clusters and GO.ID
    Clusters<-base::lapply(clusters,function(x){

      ###################
      # load the target object
      object=base::get(x[[1]])

      #################
      # check class
      if(!methods::is(object,"GO_clusters")) base::stop("object must come from ViSEAGO::GOterms_heatmap()")

      ###################
      # extract GO.ID and correspond cluster assignation to a list
      object<- methods::slot(
        methods::slot(
          object,
          "enrich_GOs"
        ),
        "data"
      )[,.(GO.ID,GO.cluster)]

      ###################
      # ordering by GO.ID
      data.table::setorder(object,GO.ID)

      ###################
      # return the object
      object$GO.cluster
    })

    ###################
    # add names
    base::names(Clusters)<-base::names(clusters)

  ###################
  # extract clusters and GO.ID
  ###################

    ###################
    # number of list elements
    n_list <- base::length(Clusters)

    ###################
    # build o matrix of one
    the_cor <- base::matrix(1, n_list, n_list)

    ###################
    # build pairwise combinations
    pairwise_combn <- utils::combn(n_list, 2)

    ###################
    # for each combination
    for (i in base::seq_len(base::ncol(pairwise_combn))) {

      ###################
      # extract l1 elements
      l1 <- pairwise_combn[1, i]

      ###################
      # extract l2 elements
      l2 <- pairwise_combn[2, i]

      ###################
      # add the value to the_cor
      the_cor[l1, l2] <- the_cor[l2, l1] <-igraph::compare(Clusters[[l1]], Clusters[[l2]],method=method)
    }

    ###################
    # add headers
    base::rownames(the_cor) <- base::colnames(the_cor) <- base::names(Clusters)

    ###################
    # return the results
    the_cor
})
