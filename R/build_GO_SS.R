#' @title build GO Semantic Similarity object.
#' @description
#' Compute the Information content (IC) on the given ontology, and
#' create a \code{\link{GO_SS-class}} object required by \code{\link{compute_SS_distances}} method to compute GO semantic similarity
#' between enriched GO terms or groups of terms.
#' @importFrom GO.db GO.db GOMFOFFSPRING GOBPOFFSPRING GOCCOFFSPRING
#' @importFrom topGO inverseList
#' @importFrom AnnotationDbi select keys as.list
#' @family GO_semantic_similarity
#' @param gene2GO a \code{\link{gene2GO-class}} object from \code{\link{annotate}} method.
#' @param enrich_GO_terms a \code{\link{enrich_GO_terms-class}} from \code{\link{merge_enrich_terms}} method.
#' @details This method use \code{\link{annotate}} and \code{\link{merge_enrich_terms}} output objects (see Arguments),
#' and compute the Information content (IC) using the internal code of \code{\link[GOSemSim]{godata}} function from \pkg{GOSemSim} package.
#' @return a \code{\link{GO_SS-class}} object required by \code{\link{compute_SS_distances}}.
#' @references
#' Alexa A, Rahnenfuhrer J, Lengauer T. Improved scoring of functional groups from gene expression data by
#' decorrelating GO graph structure. Bioinformatics 2006; 22:1600-1607.
#'
#' Guangchuang Yu, Fei Li, Yide Qin, Xiaochen Bo, Yibo Wu and Shengqi Wang. GOSemSim: an R package for measuring semantic similarity
#' among GO terms and gene products. Bioinformatics 2010 26(7):976-978.
#'
#' Herve Pages, Marc Carlson, Seth Falcon and Nianhua Li (2017). AnnotationDbi: Annotation Database Interface. R package version 1.38.0.
#' @include gene2GO.R
#' @include enrich_GO_terms.R
#' @include GO_SS.R
#' @examples
#' \dontrun{
#' # initialyse object for compute GO Semantic Similarity
#' myGOs<-ViSEAGO::build_GO_SS(
#'     myGENE2GO,
#'     BP_sResults
#' )
#' }
#' # load data example
#' utils::data(
#'     myGOs,
#'     package="ViSEAGO"
#' )
#' @name build_GO_SS
#' @rdname build_GO_SS-methods
#' @exportMethod build_GO_SS
setGeneric(
    name="build_GO_SS",
    def=function(gene2GO,enrich_GO_terms){
        standardGeneric("build_GO_SS")
})

#' @rdname build_GO_SS-methods
#' @aliases build_GO_SS
setMethod(
    "build_GO_SS",
    signature(
        gene2GO="gene2GO",
        enrich_GO_terms="enrich_GO_terms"
    ),
    definition=function(gene2GO,enrich_GO_terms){

        ## check object class

        if(!is(gene2GO,"gene2GO")){
            stop("object must be a gene2GO class object from ViSEAGO::annotate()")
        }
        if(!is(enrich_GO_terms,"enrich_GO_terms")){
            stop("object must be a enrich_GO_terms class object from ViSEAGO::merge_enrich_terms()")
        }

        ## load data

        # onto match argument
        ont<-slot(enrich_GO_terms,"ont")

        # from GO database
        go<-select(
            GO.db,
            keys=keys(GO.db),
        columns=c("GOID","ONTOLOGY")
    )

    # select a given category (MF,BP,CC)
    goids<-go$GOID[go$ONTOLOGY==ont]

    # count number of genes by term
    gocount<-lengths(
        inverseList(
            slot(gene2GO,ont)
        )
    )

    ## From GOSemsim compute IC content

    # extract gcount term names
    goname<-names(gocount)

    # ensure goterms not appearing in the specific annotation have 0 frequency
    go.diff<-setdiff(goids, goname)
    m<- double(length(go.diff))
    names(m)<-go.diff
    gocount<- as.vector(gocount)
    names(gocount)<- goname
    gocount<-c(gocount, m)

    # select offspings
    Offsprings <- switch(
        ont,
        MF =as.list(GOMFOFFSPRING),
        BP =as.list(GOBPOFFSPRING),
        CC =as.list(GOCCOFFSPRING)
    )

    # add for each term offsprings their counted values
    cnt <- gocount[goids] + vapply(goids, function(i){
        sum(gocount[Offsprings[[i]]], na.rm=TRUE)
    },0)
    names(cnt)<-goids

    # the probabilities of occurrence of GO terms in a specific corpus.
    IC<- -log(cnt/sum(gocount))

    ## return value et S4 object with compatibility with GOSemSim package

    # create object
    new(
        "GO_SS",
        db=slot(gene2GO,"db"),
        stamp =slot(gene2GO,"stamp"),
        organism=slot(gene2GO,"organism"),
        ont=ont,
        topGO=slot(enrich_GO_terms,"topGO"),
        enrich_GOs=enrich_GO_terms,
        IC = IC
    )
})
