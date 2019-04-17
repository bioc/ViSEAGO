#' @title Download EntrezGene orthologs groups.
#' @description  Download EntrezGene orthologs groups.
#' @importFrom data.table data.table
#' @details Internal function used by \code{\link{annotate}} in order to download orthologs_groups from NCBI Annotation pipeline
#'  stored in the \href{ftp://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz}{gene_group.gz} file.
#' @return a \code{\link[data.table]{data.table}}.
#' @examples
#' \dontrun{
#' ###################
#' # Organism taxid, Scientific name and common name
#' ViSEAGO::EntrezGene_orthologs()
#' }
#' @keywords internal
#' @export
EntrezGene_orthologs=function(){

  ###################
  # temp file
  temp<-base::paste(
    base::tempfile(),
    "gz",
    sep="."
  )

  ###################
  # import gene_group from NCBI
  utils::download.file(
    "ftp://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz",
    quiet=TRUE,
    destfile =temp
  )

  ##################
  #  uncompress
  R.utils::gunzip(temp)

  ##################
  # read the file
  gene_orthologs=data.table::fread(
    base::sub("\\.gz","",temp),
    verbose=FALSE,
    showProgress=FALSE
  )

  ###################
  # remove # in header
  base::names(gene_orthologs)<-base::gsub("#","",base::names(gene_orthologs))

  ###################
  # convert to character
  gene_orthologs<-gene_orthologs[,
    base::lapply(.SD,base::as.character),
    .SDcols=base::names(gene_orthologs)
  ]

  ###################
  # return the table
  base::return(gene_orthologs)
}
