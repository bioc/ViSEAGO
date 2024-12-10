#' @title perform multilevel preranked gene set enrichment analysis.
#' @description This method perform fast gene set enrichment analysis (GSEA) using \pkg{fgsea} package.
#' @import data.table
#' @import topGO
#' @import fgsea
#' @family GO_terms
#' @param geneSel a 2 columns \code{data.table} with \bold{preranked} gene identifiants (in first column) based on the statistical values (second column).
#' @param gene2GO a \code{\link{gene2GO-class}} object created by \code{\link{annotate}} method.
#' @param ont the ontology used is "MF" (Molecuar Function), "BP" (Biological Process), or "CC" (Cellular Component).
#' @param method fgsea method to use with \code{\link[fgsea]{fgseaSimple}} or \code{\link[fgsea]{fgseaMultilevel}}.
#' @param params a \code{list} with \code{\link[fgsea]{fgseaSimple}} or \code{\link[fgsea]{fgseaMultilevel}} parameters.
#' @details This method is a convenient wrapper using a given ontology category (\code{ont} argument)
#' in order to perform gene set enrichment analysis using \code{\link[fgsea]{fgseaSimple}} or
#' \code{\link[fgsea]{fgseaMultilevel}} algorithm from  \pkg{fgsea} package.
#' 
#' The complete GO annotation is required (\code{gene2GO} argument), and also a 2 columns \code{data.table} with \bold{preranked}
#' gene identifiants (in first column) based on statistical values (second column).
#' 
#' Defaults fgseaSimple parameters were used for perform test with \code{nperm} set to 10,000.\cr
#' Defaults fgseaMultilevel parameters were used for perform test 
#' except the \code{eps} arg that was set to 0 for better pvalues estimation.\cr
#' A gene frequency (\%) of leadingEdge/size is added to output \code{data.table}.
#' @return a \code{\link{fgsea-class}} object.
#' @references
#' Korotkevich G, Sukhov V, Sergushichev A (2019). "Fast gene set enrichment analysis." bioRxiv. doi: 10.1101/060012, http://biorxiv.org/content/early/2016/06/20/060012.
#' @include fgsea.R
#' @examples
#' # gene list
#' PregnantvsLactate<-data.table::fread(
#'     system.file(
#'         "extdata/data/input",
#'         "pregnantvslactate.complete.txt",
#'         package = "ViSEAGO"
#'     ),
#'     select = c("Id","padj")
#' )
#'
#' # rank Id based on statistical value (padj here)
#' PregnantvsLactate<-data.table::setorder(PregnantvsLactate,padj)
#'
#' \dontrun{
#' # connect to Bioconductor
#' Bioconductor<-ViSEAGO::Bioconductor2GO()
#'
# 'load GO annotations from Bioconductor
#' myGENE2GO<-ViSEAGO::annotate(
#'    "org.Mm.eg.db",
#'    Bioconductor
#' )
#' 
#' # run fgseaMultilevel
#' pregnantvslactate<-ViSEAGO::runfgsea(
#'     geneSel=PregnantvsLactate,
#'     gene2GO=myGENE2GO,
#'     ont="BP",
#'     method="fgseaMultilevel",
#'     params=list(
#'         minSize=5,
#'         scoreType="pos"
#'     )
#' )
#' }
#' @name runfgsea
#' @rdname runfgsea-methods
#' @exportMethod runfgsea
setGeneric(
    name="runfgsea",
    def=function(geneSel,gene2GO,ont,method=c("fgseaSimple","fgseaMultilevel"),params=list(nperm=10000,sampleSize=101,minSize=1,maxSize=Inf,eps = 0,scoreType=c("std", "pos", "neg"),nproc = 0,gseaParam = 1,BPPARAM = NULL,absEps = NULL)){
        standardGeneric("runfgsea")
    }
)

#' @rdname runfgsea-methods
#' @aliases runfgsea
setMethod(
    "runfgsea",
    signature(
        gene2GO="gene2GO",
        ont="character"
    ),
    definition=function(geneSel,gene2GO,ont,method,params){

        # fgseaSimple  default parameters
        if(method=="fgseaSimple"){
            defaults_params=list(
                nperm=10000,
                minSize=1,
                maxSize=Inf,
                scoreType=c("std", "pos", "neg"),
                nproc = 0,
                gseaParam = 1,
                BPPARAM = NULL
            )
        }

        # fgseaMultilevel default parameters
        if(method=="fgseaMultilevel"){
            defaults_params=list(
                sampleSize=101,
                minSize=1,
                maxSize=Inf,
                eps = 0,
                scoreType=c("std", "pos", "neg"),
                nproc = 0,
                gseaParam = 1,
                BPPARAM = NULL,
                absEps = NULL
            )
        }

        # params args
        params_args<-names(params)[names(params)%in%names(defaults_params)]
        
        # modify default parameters
        defaults_params[params_args]<-params[params_args]

        # rename
        params<-defaults_params

        # reverse list
        Pathways<-inverseList(slot(gene2GO,ont))

        # assign header
        names(geneSel)<-c("Id","value")

        # build ranks
        Ranks<-geneSel$value

        # assign Ids to rank
        names(Ranks)<-geneSel$Id

        # perform simple fgsea
        if(method=="fgseaSimple"){
            data<-fgseaSimple(
                Pathways,
                stats=Ranks,
                nperm=params$nperm,
                minSize=params$minSize,
                maxSize=params$maxSize,
                scoreType=params$scoreType,
                nproc=params$nproc,
                gseaParam=params$gseaParam,
                BPPARAM=params$BPPARAM
            )
        }

        # perform multilevel fgsea
        if(method=="fgseaMultilevel"){
            data<-fgseaMultilevel(
                Pathways,
                Ranks,
                sampleSize=params$sampleSize,
                minSize=params$minSize,
                maxSize=params$maxSize,
                eps=params$eps,
                scoreType=params$scoreType,
                nproc=params$nproc,
                gseaParam=params$gseaParam,
                BPPARAM=params$BPPARAM,
                absEps=params$absEps
            )
        }

        # add leadingEdge gene frequency
        data<-data[,
            `genes_frequency`:=paste(
                round(
                    length(unlist(leadingEdge))/size*100,
                    digits=3
                ),
                "% (",
                length(unlist(leadingEdge)),
                "/",
                size,
                ")",
                sep=""
            )
            ,by=pathway
        ]

        # create fgsea output
        new(
            "fgsea",
            description=paste(
                slot(gene2GO,"db"),
                slot(gene2GO,"organism"),
                slot(gene2GO,"stamp")
            ),
            ontology =ont,
            method=method,
            params=params,
            input=list(geneSel),
            data=list(data)
        )
})
