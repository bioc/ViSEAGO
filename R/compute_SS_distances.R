#' @title Compute distance between GO terms or GO clusters based on semantic similarity.
#' @description This method computes distance between GO terms or GO clusters based on semantic similarity.
#' @importFrom GO.db GO.db GOMFANCESTOR GOBPANCESTOR GOCCANCESTOR
#' @importFrom GOSemSim mgoSim
#' @importFrom AnnotationDbi select
#' @importFrom data.table data.table .I .N
#' @importFrom utils stack unstack
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
#' # load data example
#' data(
#'     myGOs,
#'     package="ViSEAGO"
#' )
#'
#' \dontrun{
#' # compute GO terms Semantic Similarity distances
#' myGOs<-ViSEAGO::compute_SS_distances(
#'     myGOs,
#'     distance=c("Resnik","Lin","Rel","Jiang","Wang")
#' )
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
#'                 minClusterSize=2
#'             )
#'         )
#'     ),
#'     samples.tree=NULL
#' )
#'
#' # compute clusters of GO terms Semantic Similarity distances
#' Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(
#'     Wang_clusters_wardD2,
#'     distance=c("max","avg","rcmax","BMA")
#' )
#' }
#' @name compute_SS_distances
#' @rdname compute_SS_distances-methods
#' @exportMethod compute_SS_distances
setGeneric(
    name="compute_SS_distances",
    def=function(object,distance){
        standardGeneric("compute_SS_distances")
    }
)

#' @rdname compute_SS_distances-methods
#' @aliases compute_SS_distances
setMethod(
    "compute_SS_distances",
    signature(
        object="ANY",
        distance="character"
    ),
    definition=function(object,distance) {

        # check object
        if(!is(object,"GO_SS") & !is(object,"GO_clusters")){

            stop("object must be a GO_SS or GO_clusters class from ViSEAGO::build_GO_SS() or ViSEAGO::GOterms_heatmap(), respectively")
        }

        if(is(object,"GO_SS")){

            # extract GO terms
            Terms<-slot(slot(object,"enrich_GOs"),"data")$GO.ID

            # for each distance
            for (x in distance){

                # check measure
                x=match.arg(x,c("Resnik","Lin","Rel","Jiang","Wang"))

                # convert GO similarity matrix to distance
                values=list(
                    as.dist(1-mgoSim(Terms,Terms,semData=object,measure=x,combine=NULL)))

                    # add distance name to values
                    names(values)=x

                    # return values
                    slot(object,"terms_dist")<-c(slot(object,"terms_dist"),values)
            }
        }else{

            # extract only selected colmuns from GOsResults
            GOclusters<-slot(slot(object,"enrich_GOs"),"data")[,c("GO.ID","term","GO.cluster"),with=FALSE]

            # cluster names
            clusters<-GOclusters[,c("GO.ID","GO.cluster"),with=FALSE]

            # convert clusters in factor
            clusters[,"GO.cluster":=factor(clusters$GO.cluster,levels = unique(clusters$GO.cluster))]

            # get ancestors
            onto=switch(
                slot(object,"ont"),
                MF=GOMFANCESTOR,
                BP=GOBPANCESTOR,
                CC=GOCCANCESTOR
            )

            # convert in list
            onto=AnnotationDbi::as.list(onto)

            # keep enrich terms
            onto=onto[clusters$GO.ID]

            # keep enrich terms
            onto=stack(onto)

            # add enrich terms
            onto=rbind(
                onto,
                cbind(
                    values=clusters$GO.ID,
                    ind=clusters$GO.ID
                )
            )

            # merge onto with GOclusters
            onto=merge(
                clusters[,c("GO.ID","GO.cluster"),with=FALSE],
                onto,
                by.x="GO.ID",
                by.y="ind",
                all.y=TRUE,
                sort=FALSE
            )

            # remove GO.ID
            onto=onto[,"GO.ID":=NULL]

            # count the term occurence by cluster
            onto=onto[,.N,by=c("GO.cluster","values")]

            # select the maxima term occurence by cluster
            onto=onto[onto[, .I[N==max(N)], by="GO.cluster"]$V1]

            # add IC to onto
            onto=data.table(onto,IC=slot(object,"IC")[onto$values])

            # select the first maxima term IC by cluster
            onto=onto[onto[, .I[which.max(IC)], by="GO.cluster"]$V1]

            # extract GO terms
            GO<-data.table(
                select(
                    GO.db,
                    keys=onto$values,
                    column="TERM"
                )
            )

            # add sizes column to GOclusters
            onto=data.table(onto,GO)

            # convert clusters in factor
            clusters<-unstack(data.frame(clusters))

            # rename clusters
            names(clusters)<-paste(onto$GO.cluster,onto$GOID,onto$TERM,sep="_")

            ## compute clusters similarities

            # number of clusters
            n=length(clusters)

            # for each distance
            for (x in distance){

                # check measure
                x=match.arg(x,c("max","avg","rcmax","BMA"))

                # create default scores matrix
                scores <- matrix(NA, nrow = n, ncol = n)

                # create default scores matrix
                for (i in seq_along(clusters)){

                    # extract the first cluster GO terms
                    GO1 <-clusters[[i]]

                    # extract second cluster GO terms
                    for (j in seq_len(i)){

                        # extract the second cluster GO terms
                        GO2<-clusters[[j]]

                        # calculate
                        scores[i, j] <-mgoSim(
                            GO1,
                            GO2,
                            semData = object,
                            measure =names(slot(object,"terms_dist")),
                            combine =x
                        )

                        # add to upper part
                        if (i != j) scores[j, i] <- scores[i, j]
                    }
                }

                # add clusters name
                colnames(scores)<-names(clusters)
                row.names(scores)<-names(clusters)

                # return matrix
                values=list(stats::as.dist(1-scores))

                # add distance name to values
                names(values)=x

                # return values
                slot(object,"clusters_dist")<-c(slot(object,"clusters_dist"),values)
            }
        }

        # return the object
        return(object)
    }
)
