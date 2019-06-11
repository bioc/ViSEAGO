#' @title Check available organisms databases at Bioconductor.
#' @description Retrieve the Bioconductor \href{http://bioconductor.org/packages/release/BiocViews.html#___OrgDb}{OrgDb}
#' available organisms databases packages.
#' @importFrom data.table data.table like
#' @importFrom AnnotationForge loadAnnDbPkgIndex
#' @importFrom AnnotationDbi species
#' @importFrom methods new
#' @family genomic_ressource
#' @details
#' This function gives genome wide annotation for available organisms databases packages from
#' \href{http://bioconductor.org/packages/release/BiocViews.html#___OrgDb}{Bioconductor OrgDb}.
#' It uses \code{loadAnnDbPkgIndex} from \pkg{AnnotationForge} package.
#' @references
#' Carlson M and Pages H (2017). AnnotationForge: Code for Building Annotation Database Packages. R package version 1.18.0.
#' @return a  \code{\link{genomic_ressource-class}} object required by \code{\link{annotate}} method.
#' @include genomic_ressource.R
#' @examples
#' ###################
#' # Check Bioconductor OrgDb available organisms
#' Bioconductor<-ViSEAGO::Bioconductor2GO()
#' @export
Bioconductor2GO=function(){

    # Available databases
    Orgs.db<-data.table(loadAnnDbPkgIndex())

    # select organisms lines and only selected columns
    Orgs.db<-Orgs.db[like(Orgs.db$Package,"^org"),c("Package","Version","species"),with=FALSE]

    # return data
    new(
        "genomic_ressource",
        db="Bioconductor",
        stamp=as.character(Sys.time()),
        data=data.table(),
        organisms=Orgs.db
    )
}
