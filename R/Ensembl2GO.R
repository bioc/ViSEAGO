#' @title Check available organisms datasets at Ensembl.
#' @description  List Ensembl referenced organisms datasets from the current (NULL) or archive (number in \code{character}) annotation version.
#' @importFrom biomaRt useEnsembl listEnsembl listDatasets
#' @importFrom data.table data.table
#' @param biomart the biomart name (eg. "ensembl", the default) available with \pkg{biomaRt} package \code{\link[biomaRt]{listEnsembl}}.
#' @param host the Ensembl host adress for \href{http://www.ensembl.org/index.html}{vertebrate} ("www.ensembl.org", the default value),
#' \href{http://www.plants.ensembl.org/index.html}{plants} ("plants.ensembl.org"),
#' \href{http://www.metazoa.ensembl.org/index.html}{metazoa} ("metazoa.ensembl.org"),
#'  or \href{http://www.fungi.ensembl.org/index.html}{fungi} ("fungi.ensembl.org").
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
#' # host address to use for Ensembl vertebrate/Plants/Protists/Metazoa/Fungi
#' # host="www.ensembl.org" # vertebrate
#' # host="plants.ensembl.org"
#' # host="protists.ensembl.org"
#' # host="metazoa.ensembl.org"
#' # host="fungi.ensembl.org"
#' # biomart is not available for bacteria
#'
#' \dontrun{
#' # check the Ensembl available biomart (if not known)
#' # for Animals (host="www.ensembl.org", the default)
#'  biomaRt::listEnsembl()
#'
#' # List Ensembl available organisms
#' Ensembl<-ViSEAGO::Ensembl2GO(
#'  biomart="ensembl",
#'  host="www.ensembl.org",
#'  version=NULL
#' )
#' }
#' @export
Ensembl2GO=function(biomart="ensembl",host="www.ensembl.org",version=NULL){

    # check the ensembl genes releases
    Ensembl<-listEnsembl(host=host,version=version)

    # check the ensembl host versus mart name
    match.arg(biomart,Ensembl$biomart)

    # connect to Ensembl
    mart<-useEnsembl(biomart,host=host,version=version)

    # return data in genomic_ressource class
    new(
        "genomic_ressource",
        db="Ensembl",
        stamp=paste(host,Ensembl$version[Ensembl$biomart==biomart]),
        data=data.table(),
        mart=list(mart),
        organisms=data.table(listDatasets(mart))
    )
}
