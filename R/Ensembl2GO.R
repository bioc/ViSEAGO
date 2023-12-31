#' @title Check available organisms datasets at Ensembl.
#' @description  List Ensembl referenced organisms datasets from the current (NULL) or archive (number in \code{character}) annotation version.
#' @importFrom biomaRt useEnsembl listEnsembl listEnsemblGenomes useEnsemblGenomes listDatasets 
#' @importFrom data.table data.table
#' @param biomart the biomart name available with \pkg{biomaRt} package \code{\link[biomaRt]{listEnsembl}} ("genes", the default) or \code{\link[biomaRt]{listEnsemblGenomes}} ("protists_mart", "fungi_mart", "plants_mart").
#' @param GRCh GRCh version to connect to if not the current GRCh38, currently this can only be 37
#' @param version the annotation version to use (eg. NULL for the default current version, or a version number in \code{character})
#' @family genomic_ressource
#' @details
#' This function gives referenced organisms genomes at \href{http://ensemblgenomes.org/}{Ensembl}.
#' It uses the \code{\link[biomaRt]{useEnsembl}} and  \code{\link[biomaRt]{listDatasets}} from \pkg{biomaRt} package.
#' @return a  \code{\link{genomic_ressource-class}} object required by \code{\link{annotate}}.
#' @references
#' Durinck S, Spellman P, Birney E and Huber W (2009).
#' Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt.
#' Nature Protocols, 4, pp. 1184-1191.
#'
#' Durinck S, Moreau Y, Kasprzyk A, Davis S, De Moor B, Brazma A and Huber W (2005).
#' BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis. Bioinformatics, 21, pp. 3439-3440.
#'
#' Matt Dowle and Arun Srinivasan (2017). data.table: Extension of data.frame. R package version 1.10.4. https://CRAN.R-project.org/package=data.table.
#' @include genomic_ressource.R
#' @examples
#' \dontrun{
#' # check the Ensembl available biomart (if not known)
#' biomaRt::listEnsembl()
#'
#' # List Ensembl available organisms
#' Ensembl<-ViSEAGO::Ensembl2GO(
#'  biomart="genes",
#'  GRCh = NULL,
#'  version=NULL
#' )
#' }
#' @export
Ensembl2GO=function(biomart="genes",GRCh = NULL,version=NULL){

    # check the ensembl host versus mart name
    match.arg(biomart,c("genes","protists_mart","fungi_mart","metazoa_mart","plants_mart"))
    
    # Connect to Ensembl
    if(biomart=="genes"){

        # check the ensembl genes releases
        Ensembl<-listEnsembl(
            GRCh =GRCh,
            version=version
        )

        # connect to Ensembl
        mart<-useEnsembl(
            biomart,
            GRCh =GRCh,
            version=version
        )

    }else{

        # check the ensembl genes releases
        Ensembl<-listEnsemblGenomes()

        # connect to Ensembl
        mart<-useEnsemblGenomes(biomart)
    }

    # return data in genomic_ressource class
    new(
        "genomic_ressource",
        db="Ensembl",
        stamp=paste(biomart,GRCh,Ensembl$version[Ensembl$biomart==biomart]),
        data=data.table(),
        mart=list(mart),
        organisms=data.table(listDatasets(mart))
    )
}
