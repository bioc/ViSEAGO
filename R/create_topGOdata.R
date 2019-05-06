#' @title Create topGOdata object for enrichment test with topGO package.
#' @description This method create a \code{\link[topGO]{topGOdata-class}} object required by \pkg{topGO} package in order
#' to perform GO enrichment test.
#' @import topGO
#' @family GO_terms
#' @param geneSel genes of interest.
#' @param allGenes customized background genes.
#' @param geneList logical factor (1: genes of interest, 0: genes background,
#' and gene identifiants in \code{names}) (default value to NULL).
#' @param gene2GO a \code{\link{gene2GO-class}} object created by \code{\link{annotate}} method.
#' @param ont the ontology used is "MF" (Molecuar Function), "BP" (Biological Process), or "CC" (Cellular Component).
#' @param nodeSize the minimum number of genes for each GO term.
#' @details This method is a convenient wrapper building a \code{\link[topGO]{topGOdata-class}} object
#' using a given ontology category (\code{ont} argument) in order to perform GO enrichment test.
#' The complete GO annotation is required (\code{gene2GO} argument) and also
#' the list of genes of interest (\code{geneSel} argument) against the corresponding background (\code{allGenes} argument)
#' separately, or grouped together in a \code{factor} (\code{geneList} argument).
#' @return a  \code{\link[topGO]{topGOdata-class}} object required by \code{runTest} from \pkg{topGO} package.
#' @references
#' Alexa A, Rahnenfuhrer J, Lengauer T. Improved scoring of functional groups from gene expression data by
#' decorrelating GO graph structure. Bioinformatics 2006; 22:1600-1607.
#' @include gene2GO.R
#' @examples
#'  ###################
#'  # load genes identifiants (GeneID,ENS...) background (Expressed genes)
#'  background<-scan(
#'   system.file(
#'    "extdata/data/input",
#'    "background_L.txt",
#'    package = "ViSEAGO"
#'   ),
#'   quiet=TRUE,
#'   what=""
#'  )
#'
#'  ###################
#'  # load Differentialy Expressed (DE) gene identifiants from files
#'  pregnantvslactateDE<-scan(
#'   system.file(
#'    "extdata/data/input",
#'    "pregnantvslactateDE.txt",
#'    package = "ViSEAGO"
#'  ),
#'   quiet=TRUE,
#'   what=""
#'  )
#' \dontrun{
#' ###################
#' # create topGOdata for BP for each list of DE genes
#' BP_L_pregnantvslactate<-ViSEAGO::create_topGOdata(
#'  geneSel=pregnantvslactateDE,
#'  allGenes=background,
#'  gene2GO=myGENE2GO,
#'  ont="BP",
#'  nodeSize=5
#' )
#' }
#' @name create_topGOdata
#' @rdname create_topGOdata-methods
#' @exportMethod create_topGOdata
setGeneric(
    name="create_topGOdata",
    def=function(geneSel,allGenes,geneList=NULL,gene2GO,ont,nodeSize){
        standardGeneric("create_topGOdata")
    }
)

#' @rdname create_topGOdata-methods
#' @aliases create_topGOdata
setMethod(
    "create_topGOdata",
    signature(
        gene2GO="gene2GO",
        ont="character",
        nodeSize="numeric"
    ),
    definition=function(geneSel,allGenes,geneList,gene2GO,ont,nodeSize){

        ## check list


        # check object class
        if(!is(gene2GO,"gene2GO")){
            stop("object must be a gene2GO class obtained using ViSEAGO::annotate()")
        }

        # check ontology
        ont=match.arg(ont,c("MF","BP","CC"))

        ## build allGenes factor

        if(is.null(geneList)){

            # check duplicate in allGenes
            if(length(allGenes)!=length(unique(allGenes))){

                # warning message
                message("allGenes contain genes redondancy.\nduplicate elements were removed.")

                # remove duplicate
                allGenes=unique(allGenes)
            }

            # check duplicate in geneSel
            if(length(geneSel)!=length(unique(geneSel))){

                # warning message
                message("geneSel contain genes redondancy.\nduplicate elements were removed.")

                # remove duplicate
                geneSel=unique(geneSel)
            }

            # create geneList object
            AllGenes<-rep(0,length(allGenes))

            # assign genes_ref names to geneList
            attr(AllGenes,"names")<-allGenes

            # assign value 1 for DE genes in geneList
            AllGenes[attr(AllGenes,"names")%in%geneSel]<-1

            # convert geneList in factor (used by topGO)
            AllGenes<-factor(AllGenes)

        }else{

            # check duplicate in geneSel
            if(length(names(geneList))!=length(unique(names(geneList)))){

                # warning message
                stop("geneList contain genes redondancy.")
            }

            # rename geneList
            AllGenes<-geneList
        }

        # load topGO
        requireNamespace("topGO")

        # create GOdata
        new(
            "topGOdata",
            description=paste(
                slot(gene2GO,"db"),
                slot(gene2GO,"organism"),
                slot(gene2GO,"stamp")
            ),
            ontology =ont,
            allGenes = AllGenes,
            annot = annFUN.gene2GO,
            nodeSize =nodeSize,
            gene2GO =slot(gene2GO,ont)
        )
    }
)
