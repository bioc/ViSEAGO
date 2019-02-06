#' @title Compute distance between GO terms or GO clusters based on semantic similarity.
#' @description This method computes distance between GO terms or GO clusters based on semantic similarity.
#' @importFrom methods setGeneric setMethod
#' @importFrom GO.db GOMFANCESTOR GOBPANCESTOR GOCCANCESTOR
#' @importFrom GOSemSim mgoSim
#' @importFrom AnnotationDbi select
#' @importFrom data.table data.table
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
#' # compute GO terms Semantic Similarity distances
#' myGOs<-ViSEAGO::compute_SS_distances(
#'  myGOs,
#'  distance=c("Resnik","Lin","Rel","Jiang","Wang")
#' )
#'
#' ###################
#' # compute clusters of GO terms Semantic Similarity distances
#' Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(
#'  Wang_clusters_wardD2,
#'  distance=c("max","avg","rcmax","BMA")
#' )
#' @export
setGeneric(name="compute_SS_distances",def=function(object,distance) {standardGeneric("compute_SS_distances")})

setMethod("compute_SS_distances",definition=function(object,distance) {

  ##################
  # check object
  if(!base::class(object)%in%c("GO_SS","GO_clusters"))base::stop("object must be a GO_SS or GO_clusters class from ViSEAGO::build_GO_SS() or ViSEAGO::GOterms_heatmap(), respectively")

  if(base::class(object)=="GO_SS"){

    ###################
    # extract GO terms
    Terms<-object@enrich_GOs@data$GO.ID

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
      object@terms_dist<-base::c(object@terms_dist,values)
    }
  }else{

    ###################
    # extract only selected colmuns from GOsResults
    GOclusters<-object@enrich_GOs@data[,.(GO.ID,term,GO.cluster)]

    ###################
    # cluster names
    clusters<-GOclusters[,.(GO.ID,GO.cluster)]

    ###################
    # convert clusters in factor
    clusters[,GO.cluster:=factor(GO.cluster,levels = unique(GO.cluster))]

    ###################
    # get ancestors
    onto=base::switch(object@ont,MF=GO.db::GOMFANCESTOR,BP=GO.db::GOBPANCESTOR,CC=GO.db::GOCCANCESTOR)

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
    onto=base::rbind(onto,base::cbind(values=clusters$GO.ID,ind=clusters$GO.ID))

    ###################
    # merge onto with GOclusters
    onto=merge(clusters[,.(GO.ID,GO.cluster)],onto,by.x="GO.ID",by.y="ind",all.y=T,sort=F)

    ###################
    # remove GO.ID
    onto=onto[,GO.ID:=NULL]

    ###################
    # count the term occurence by cluster
    onto=onto[,.N,by=list(GO.cluster,values)]

    ###################
    # select the maxima term occurence by cluster
    onto=onto[onto[, .I[N== max(N)], by=GO.cluster]$V1]

    ###################
    # add IC to onto
    onto=data.table::data.table(onto,IC=object@IC[onto$values])

    ###################
    # select the first maxima term IC by cluster
    onto=onto[onto[, .I[which.max(IC)], by=GO.cluster]$V1]

    ###################
    # extract GO terms
    GO<-data.table::data.table(AnnotationDbi::select(GO.db::GO.db, keys=onto$values,column ="TERM"))

    ###################
    #  add sizes column to GOclusters
    onto=data.table::data.table(onto,GO)

    ###################
    # convert clusters in factor
    clusters<-utils::unstack(data.frame(clusters))

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
        for (i in base::seq_along(clusters)) {

          ###################
          # extract the first cluster GO terms
          GO1 <-clusters[[i]]

          ###################
          # extract second cluster GO terms
          for (j in 1:i) {

            ###################
            # extract the second cluster GO terms
            GO2<-clusters[[j]]

            ###################
            # calculate
            scores[i, j] <-GOSemSim::mgoSim(GO1,GO2,semData = object,
            measure =base::names(object@terms_dist), combine =x)

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
        object@clusters_dist<-c(object@clusters_dist,values)
    }
  }

  ###################
  # return the object
  return(object)
})
