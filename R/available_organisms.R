#' @title Display available organisms from a specified database.
#' @description Display an interactive table with available organisms from a genomic ressource database (Bioconductor, EntrezGene, Ensembl, Uniprot).
#' @importFrom DT datatable
#' @family genomic_ressource
#' @family visualization
#' @param object a  \code{\link{genomic_ressource-class}} object created by \code{\link{Bioconductor2GO}}, \code{\link{EntrezGene2GO}},
#' \code{\link{Ensembl2GO}},or \code{\link{Uniprot2GO}} methods.
#' @details
#' an interactive \code{datatable}.
#' @return javascript datatable
#' @references
#'  Yihui Xie (2016). DT: A Wrapper of the JavaScript Library 'DataTables'. R package version 0.2. https://CRAN.R-project.org/package=DT
#' @include genomic_ressource.R
#' @examples
#' ###################
#' # display Bioconductor table
#' Bioconductor<-ViSEAGO::Bioconductor2GO()
#' ViSEAGO::available_organisms(Bioconductor)
#'
#' ###################
#' # display EntrezGene table
#' EntrezGene<-ViSEAGO::EntrezGene2GO()
#' ViSEAGO::available_organisms(EntrezGene)
#'
#' ###################
#' # display Ensembl table
#' Ensembl<-ViSEAGO::Ensembl2GO()
#' ViSEAGO::available_organisms(Ensembl)
#'
#' ###################
#' # display Uniprot table
#' Uniprot<-ViSEAGO::Uniprot2GO()
#' ViSEAGO::available_organisms(Uniprot)
#'
#' @export
setGeneric(name="available_organisms",def=function(object) {standardGeneric("available_organisms")})

setMethod("available_organisms",signature="genomic_ressource",definition=function(object) {

  ##################
  # create a datatable
  DT::datatable(object@organisms,

    ###################
    # table width
    width=650,

    ###################
    # filter column
    filter ='top',

    ###################
    # extensions to use
    extensions =c('Scroller'),

    ###################
    # DT options to use
    options = base::list(

      ###################
      # the dom
      dom='ltipr',

      ###################
      # page length
      pageLength =5,

      ###################
      # y scroller
      deferRender = TRUE,
      scrollY = 200,
      scroller = TRUE
    )
  )
})
