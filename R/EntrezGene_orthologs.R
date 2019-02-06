#' @title Download EntrezGene orthologs groups.
#' @description  Download EntrezGene orthologs groups.
#' @importFrom data.table data.table
#' @details Internal function used by \code{\link{annotate}} in order to download orthologs_groups from NCBI Annotation pipeline
#'  stored in the \href{ftp://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz}{gene_group.gz} file.
#' @return a \code{\link[data.table]{data.table}}.
#' @examples
#' ###################
#' # Organism taxid, Scientific name and common name
#' ViSEAGO::EntrezGene_orthologs()
#' @keywords internal
#' @export
EntrezGene_orthologs=function(){

  ###################
  # import gene_group from NCBI
  utils::download.file("ftp://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz",quiet=T,
  destfile = "gene_orthologs.gz")

  ##################
  #  uncompress
  R.utils::gunzip("gene_orthologs.gz")

  ##################
  # read the file
  gene_orthologs=data.table::fread("gene_orthologs",verbose=F,showProgress=F)

  ###################
  # unlink gz file
  base::unlink("gene_orthologs.gz")
  base::unlink("gene_orthologs")

  ###################
  # remove # in header
  base::names(gene_orthologs)<-base::gsub("#","",base::names(gene_orthologs))

  ###################
  # return the table
  base::return(gene_orthologs)
}
