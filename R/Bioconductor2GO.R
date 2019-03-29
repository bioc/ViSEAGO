#' @title Check available organisms databases at Bioconductor.
#' @description Retrieve the Bioconductor \href{http://bioconductor.org/packages/release/BiocViews.html#___OrgDb}{OrgDb}
#' available organisms databases packages.
#' @importFrom data.table data.table like
#' @importFrom methods new
#' @importFrom AnnotationForge loadAnnDbPkgIndex
#' @importFrom AnnotationDbi species
#' @family genomic_ressource
#' @details
#' This function gives genome wide annotation for available organisms databases packages from
#' \href{http://bioconductor.org/packages/release/BiocViews.html#___OrgDb}{Bioconductor OrgDb}.
#' It uses \code{\link[AnnotationForge]{loadAnnDbPkgIndex}} from \pkg{AnnotationForge} package.
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

  ###################
  # Available databases
  Orgs.db<-data.table::data.table(AnnotationForge::loadAnnDbPkgIndex())

  ###################
  # select organisms lines and only selected columns
  Orgs.db<-Orgs.db[data.table::like(Package,"^org"),.(Package,Version,species)]

  ###################
  # return data
  methods::new("genomic_ressource",
               db="Bioconductor",
               stamp=base::as.character(base::Sys.time()),
               data=data.table::data.table(),
               organisms=Orgs.db
              )
}
