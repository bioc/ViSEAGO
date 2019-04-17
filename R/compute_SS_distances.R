#' @title Compute distance between GO terms or GO clusters based on semantic similarity.
#' @description This method computes distance between GO terms or GO clusters based on semantic similarity.
#' @importFrom methods setGeneric setMethod slot is signature
#' @importFrom GO.db GOMFANCESTOR GOBPANCESTOR GOCCANCESTOR
#' @importFrom GOSemSim mgoSim
#' @importFrom AnnotationDbi select
#' @importFrom data.table data.table .I
#' @family  GO_semantic_similarity
#' @param object a \code{\link{GO_SS-class}}, or  \code{\link{GO_clusters-class}} objects created by \code{\link{build_GO_SS}}
#'  or \code{\link{GOterms_heatmap}} methods, respectively.
#' @param distance The available methods for calculating GO terms Semantic Similarity (SS) are
#' "Resnik", "Rel", "Lin", and "Jiang" which are based on Information Content (IC), and "Wang" which is based on graph topology.\cr
#' The available methods for calculating clusters of GO terms SS are "max", "avg","rcmax", and "BMA".
#' @details This method computes semantic similarity distances between all GO terms
#' provided by \code{\link{GO_SS-class}} object.\cr
#' This method also computes semantic similarity distances between all GO clusters
#' provided by \code{\link{GO_clusters-class}} object.\cr
#'
#' Semantic Similarity computations are based on \code{\link[GOSemSim]{mgoSim}} method from the \pkg{GoSemSim} package.
#' @return a \code{\link{GO_SS-class}}, or a \code{\link{GO_clusters-class}} object (same class as input object).
#' @references
#'
#' Marc Carlson (2017). GO.db: A set of annotation maps describing the entire Gene Ontology. R package version 3.4.1.
#'
#' Guangchuang Yu, Fei Li, Yide Qin, Xiaochen Bo, Yibo Wu and Shengqi Wang. GOSemSim: an R package for measuring semantic similarity
#' among GO terms and gene products. Bioinformatics 2010 26(7):976-978
#'
#' Herve Pages, Marc Carlson, Seth Falcon and Nianhua Li (2017). AnnotationDbi: Annotation Database Interface. R package version 1.38.0.
#' @include GO_SS.R GO_clusters.R
#' @examples
#' ###################
#' # load data example
#' utils::data(
#'  myGOs,
#'  package="ViSEAGO"
#' )
#'
#' \dontrun{
#' ###################
#' # compute GO terms Semantic Similarity distances
#' myGOs<-ViSEAGO::compute_SS_distances(
#'  myGOs,
#'  distance=c("Resnik","Lin","Rel","Jiang","Wang")
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
#'  distance=c("max","avg","rcmax","BMA")
#' )
#' }
#' @name compute_SS_distances
#' @rdname compute_SS_distances-methods
#' @exportMethod compute_SS_distances
setGeneric(name="compute_SS_distances",def=function(object,distance){standardGeneric("compute_SS_distances")})

#' @rdname compute_SS_distances-methods
#' @aliases compute_SS_distances
setMethod("compute_SS_distances",
  methods::signature(
    object="ANY",
    distance="character"
  ),definition=function(object,distance) {

  ##################
  # check object
  if(!base::class(object)%in%c("GO_SS","GO_clusters")){
    base::stop("object must be a GO_SS or GO_clusters class from ViSEAGO::build_GO_SS() or ViSEAGO::GOterms_heatmap(), respectively")
  }

  if(methods::is(object,"GO_SS")){

    ###################
    # extract GO terms
    Terms<-methods::slot(
      methods::slot(
        object,
        "enrich_GOs"
      ),
      "data"
    )$GO.ID

    ###################
    # for each distance
    for (x in distance){

      ###################
      # check measure
      x=base::match.arg(x,c("Resnik","Lin","Rel","Jiang","Wang"))

      ###################
      # convert GO similarity matrix to distance
      values=base::list(as.dist(1-GOSemSim::mgoSim(
      Terms,Terms,semData=object,measure=x,combine=NULL)))

      ###################
      # add distance name to values
      names(values)=x

      ###################
      # return values
      methods::slot(object,"terms_dist")<-base::c(methods::slot(object,"terms_dist"),values)
    }
  }else{

    ###################
    # extract only selected colmuns from GOsResults
    GOclusters<-methods::slot(
      methods::slot(
        object,
        "enrich_GOs"
      ),
      "data"
    )[,base::c("GO.ID","term","GO.cluster"),with=FALSE]

    ###################
    # cluster names
    clusters<-GOclusters[,base::c("GO.ID","GO.cluster"),with=FALSE]

    ###################
    # convert clusters in factor
    clusters[,"GO.cluster":=base::factor(
      clusters$GO.cluster,
      levels = base::unique(clusters$GO.cluster)
    )]

    ###################
    # get ancestors
    onto=base::switch(
      methods::slot(object,"ont"),
      MF=GO.db::GOMFANCESTOR,
      BP=GO.db::GOBPANCESTOR,
      CC=GO.db::GOCCANCESTOR
    )

    ###################
    # convert in list
    onto=AnnotationDbi::as.list(onto)

    ###################
    # keep enrich terms
    onto=onto[clusters$GO.ID]

    ###################
    # keep enrich terms
    onto=utils::stack(onto)

    ###################
    # add enrich terms
    onto=base::rbind(
      onto,
      base::cbind(
        values=clusters$GO.ID,
        ind=clusters$GO.ID
      )
    )

    ###################
    # merge onto with GOclusters
    onto=merge(
      clusters[,base::c("GO.ID","GO.cluster"),with=FALSE],
      onto,
      by.x="GO.ID",
      by.y="ind",
      all.y=T,
      sort=F
    )

    ###################
    # remove GO.ID
    onto=onto[,"GO.ID":=NULL]

    ###################
    # count the term occurence by cluster
    onto=onto[,.N,by=c("GO.cluster","values")]

    ###################
    # select the maxima term occurence by cluster
    onto=onto[onto[, .I[N==base::max(N)], by="GO.cluster"]$V1]

    ###################
    # add IC to onto
    onto=data.table::data.table(onto,IC=object@IC[onto$values])

    ###################
    # select the first maxima term IC by cluster
    onto=onto[onto[, .I[base::which.max(IC)], by="GO.cluster"]$V1]

    ###################
    # extract GO terms
    GO<-data.table::data.table(
      AnnotationDbi::select(
        GO.db::GO.db,
        keys=onto$values,
        column ="TERM"
      )
    )

    ###################
    #  add sizes column to GOclusters
    onto=data.table::data.table(onto,GO)

    ###################
    # convert clusters in factor
    clusters<-utils::unstack(base::data.frame(clusters))

    ###################
    # rename clusters
    base::names(clusters)<-paste(onto$GO.cluster,onto$GOID,onto$TERM,sep="_")

    ###################
    # compute clusters similarities
    ###################

      ###################
      # number of clusters
      n=base::length(clusters)

      ###################
      # for each distance
      for (x in distance){

        ###################
        # check measure
        x=base::match.arg(x,base::c("max","avg","rcmax","BMA"))

        ###################
        # create default scores matrix
        scores <- base::matrix(NA, nrow = n, ncol = n)

        ###################
        # create default scores matrix
        for (i in base::seq_along(clusters)){

          ###################
          # extract the first cluster GO terms
          GO1 <-clusters[[i]]

          ###################
          # extract second cluster GO terms
          for (j in base::seq_len(i)){

            ###################
            # extract the second cluster GO terms
            GO2<-clusters[[j]]

            ###################
            # calculate
            scores[i, j] <-GOSemSim::mgoSim(
              GO1,
              GO2,
              semData = object,
              measure =base::names(methods::slot(object,"terms_dist")),
              combine =x
            )

            ###################
            # add to upper part
            if (i != j) scores[j, i] <- scores[i, j]
          }
        }

        ###################
        # add clusters name
        base::colnames(scores)<-base::names(clusters)
        base::row.names(scores)<-base::names(clusters)

        ###################
        # return matrix
        values=base::list(stats::as.dist(1-scores))

        ###################
        # add distance name to values
        base::names(values)=x

        ###################
        # return values
        methods::slot(object,"clusters_dist")<-c(methods::slot(object,"clusters_dist"),values)
    }
  }

  ###################
  # return the object
  return(object)
})
