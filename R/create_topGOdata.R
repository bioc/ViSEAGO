#' @title Create topGOdata object for enrichment test with topGO package.
#' @description This method create a \code{\link[topGO]{topGOdata-class}} object required by \pkg{topGO} package in order
#' to perform GO enrichment test.
#' @importFrom methods setGeneric setMethod new slot is
#' @import topGO
#' @importFrom topGO annFUN.gene2GO
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
#' @return a  \code{\link[topGO]{topGOdata-class}} object required by \code{\link[topGO]{runTest}}.
#' @references
#' Alexa A, Rahnenfuhrer J, Lengauer T. Improved scoring of functional groups from gene expression data by
#' decorrelating GO graph structure. Bioinformatics 2006; 22:1600-1607.
#' @examples
#' \dontrun{
#'  ###################
#'  # load genes identifiants (GeneID,ENS...) background (Expressed genes)
#'  background_L<-base::scan(
#'   base::system.file("extdata/data/input","background_L.txt",package = "ViSEAGO"),
#'   quiet=T,what=""
#'  )
#'
#'  ###################
#'  # load Differentialy Expressed (DE) gene identifiants from files
#'  L_pregnantvslactateDE<-base::scan(
#'   base::system.file("extdata/data/input","L_pregnantvslactateDE.txt",package = "ViSEAGO"),
#'   quiet=T,what=""
#'  )
#'
#'  L_virginvslactateDE<-base::scan(
#'   base::system.file("extdata/data/input","L_virginvslactateDE.txt",package = "ViSEAGO"),
#'   quiet=T,what=""
#'  )
#'
#'  L_virginvspregnantDE<-base::scan(
#'   base::system.file("extdata/data/input","L_virginvspregnantDE.txt",package = "ViSEAGO"),
#'   quiet=T,what=""
#'  )
#'
#' ###################
#' # create topGOdata for BP for each list of DE genes
#' BP_L_pregnantvslactate<-ViSEAGO::create_topGOdata(
#'  geneSel=L_pregnantvslactateDE,
#'  allGenes=background_L,
#'  gene2GO=myGENE2GO,
#'  ont="BP",
#'  nodeSize=5
#' )
#'
#' BP_L_virginvslactate<-ViSEAGO::create_topGOdata(
#'  geneSel=L_virginvslactateDE,
#'  allGenes=background_L,
#'  gene2GO=myGENE2GO,
#'  ont="BP",
#'  nodeSize=5
#' )
#'
#' BP_L_virginvspregnant<-ViSEAGO::create_topGOdata(
#'  geneSel=L_virginvspregnantDE,
#'  allGenes=background_L,
#'  gene2GO=myGENE2GO,
#'  ont="BP",
#'  nodeSize=5
#' )
#'
#' ###################
#' # perform topGO tests
#' elim_BP_L_pregnantvslactate<-topGO::runTest(
#'  BP_L_pregnantvslactate,
#'  algorithm ="elim",
#'  statistic = "fisher"
#' )
#'
#' elim_BP_L_virginvslactate<-topGO::runTest(
#'  BP_L_virginvslactate,
#'  algorithm ="elim",
#'  statistic = "fisher"
#' )
#'
#' elim_BP_L_virginvspregnant<-topGO::runTest(
#'  BP_L_virginvspregnant,
#'  algorithm ="elim",
#'  statistic = "fisher"
#' )
#' }
#' @exportMethod create_topGOdata
setGeneric(name="create_topGOdata",def=function(geneSel,allGenes,geneList=NULL,gene2GO,ont,nodeSize){standardGeneric("create_topGOdata")})

setMethod("create_topGOdata",definition=function(geneSel,allGenes,geneList,gene2GO,ont,nodeSize){

  ###################
  # check list
  ###################

    ###################
    # check object class
    if(!methods::is(gene2GO,"gene2GO"))base::stop("object must be a GENE2GO class obtained using ViSEAGO::annotate()")

    ###################
    # check ontology
    ont=base::match.arg(ont,c("MF","BP","CC"))

  ###################
  # build allGenes factor
  ###################

    if(base::is.null(geneList)){

      ###################
      # check duplicate in allGenes
      if(base::length(allGenes)!=base::length(base::unique(allGenes))){

        ###################
        # warning message
        base::warning("allGenes contain genes redondancy.\nduplicate elements were removed.")

        ###################
        # remove duplicate
        allGenes=base::unique(allGenes)
      }

      ###################
      # check duplicate in geneSel
      if(base::length(geneSel)!=base::length(base::unique(geneSel))){

        ###################
        # warning message
        base::warning("geneSel contain genes redondancy.\nduplicate elements were removed.")

        ###################
        # remove duplicate
        geneSel=base::unique(geneSel)
      }

      ###################
      # create geneList object
      AllGenes<-base::rep(0,base::length(allGenes))

      ###################
      # assign genes_ref names to geneList
      base::attr(AllGenes,"names")<-allGenes

      ###################
      # assign value 1 for DE genes in geneList
      AllGenes[base::attr(AllGenes,"names")%in%geneSel]<-1

      ###################
      # convert geneList in factor (used by topGO)
      AllGenes<-base::factor(AllGenes)

    }else{

      ###################
      # check duplicate in geneSel
      if(base::length(base::names(geneList))!=base::length(base::unique(base::names(geneList)))){

        ###################
        # warning message
        base::stop("geneList contain genes redondancy.")
      }

      ###################
      # rename geneList
      AllGenes<-geneList
  }

  ###################
  # topGO
  ###################

    ###################
    # load topGO package
    base::require("topGO")

    ###################
    # create GOdata
    methods::new(
      "topGOdata",
      description=base::paste(
        methods::slot(gene2GO,"db"),
        methods::slot(gene2GO,"organism"),
        methods::slot(gene2GO,"stamp")
      ),
      ontology =ont,
      allGenes = AllGenes,
      annot = topGO::annFUN.gene2GO,
      nodeSize =nodeSize,
      gene2GO =methods::slot(gene2GO,ont)
    )
})
