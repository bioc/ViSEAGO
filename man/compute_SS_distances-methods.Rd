% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_SS_distances.R
\name{compute_SS_distances}
\alias{compute_SS_distances}
\alias{compute_SS_distances,ANY,character-method}
\title{Compute distance between GO terms or GO clusters based on semantic similarity.}
\usage{
compute_SS_distances(object, distance)

\S4method{compute_SS_distances}{ANY,character}(object, distance)
}
\arguments{
\item{object}{a \code{\link{GO_SS-class}}, or  \code{\link{GO_clusters-class}} objects created by \code{\link{build_GO_SS}}
or \code{\link{GOterms_heatmap}} methods, respectively.}

\item{distance}{The available methods for calculating GO terms Semantic Similarity (SS) are
"Resnik", "Rel", "Lin", and "Jiang" which are based on Information Content (IC), and "Wang" which is based on graph topology.\cr
The available methods for calculating clusters of GO terms SS are "max", "avg","rcmax", and "BMA".}
}
\value{
a \code{\link{GO_SS-class}}, or a \code{\link{GO_clusters-class}} object (same class as input object).
}
\description{
This method computes distance between GO terms or GO clusters based on semantic similarity.
}
\details{
This method computes semantic similarity distances between all GO terms
provided by \code{\link{GO_SS-class}} object.\cr
This method also computes semantic similarity distances between all GO clusters
provided by \code{\link{GO_clusters-class}} object.\cr

Semantic Similarity computations are based on \code{\link[GOSemSim]{mgoSim}} method from the \pkg{GoSemSim} package.
}
\examples{
# load data example
data(
    myGOs,
    package="ViSEAGO"
)

\dontrun{
# compute GO terms Semantic Similarity distances
myGOs<-ViSEAGO::compute_SS_distances(
    myGOs,
    distance=c("Resnik","Lin","Rel","Jiang","Wang")
)

# GOtermsHeatmap with default parameters
Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
    myGOs,
    showIC=TRUE,
    showGOlabels=TRUE,
    GO.tree=list(
        tree=list(
            distance="Wang",
            aggreg.method="ward.D2",
            rotate=NULL
        ),
        cut=list(
            dynamic=list(
                pamStage=TRUE,
                pamRespectsDendro=TRUE,
                deepSplit=2,
                minClusterSize=2
            )
        )
    ),
    samples.tree=NULL
)

# compute clusters of GO terms Semantic Similarity distances
Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(
    Wang_clusters_wardD2,
    distance=c("max","avg","rcmax","BMA")
)
}
}
\references{
Marc Carlson (2017). GO.db: A set of annotation maps describing the entire Gene Ontology. R package version 3.4.1.

Guangchuang Yu, Fei Li, Yide Qin, Xiaochen Bo, Yibo Wu and Shengqi Wang. GOSemSim: an R package for measuring semantic similarity
among GO terms and gene products. Bioinformatics 2010 26(7):976-978

Herve Pages, Marc Carlson, Seth Falcon and Nianhua Li (2017). AnnotationDbi: Annotation Database Interface. R package version 1.38.0.
}
\seealso{
Other GO_semantic_similarity: 
\code{\link{GO_SS-class}},
\code{\link{build_GO_SS}()}
}
\concept{GO_semantic_similarity}
