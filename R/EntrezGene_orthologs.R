#' @title Download EntrezGene orthologs groups.
#' @description  Download EntrezGene orthologs groups.
#' @importFrom data.table data.table fread .SD
#' @importFrom R.utils gunzip
#' @details Internal function used by \code{\link{annotate}} in order to download orthologs_groups from NCBI Annotation pipeline
#'  stored in the \href{ftp://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz}{gene_group.gz} file.
#' @return a \code{\link[data.table]{data.table}}.
#' @examples
#' \dontrun{
#' # Organism taxid, Scientific name and common name
#' ViSEAGO::EntrezGene_orthologs()
#' }
#' @keywords internal
#' @export
EntrezGene_orthologs=function(){

    # temp file
    temp<-paste(
        tempfile(),
        "gz",
        sep="."
    )

    # import gene_group from NCBI
    download.file(
        "https://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz",
        quiet=TRUE,
        destfile =temp,
        method = "libcurl"
    )

    # read the file
    gene_orthologs=fread(
        temp,
        verbose=FALSE,
        showProgress=FALSE
    )

    # remove # in header
    names(gene_orthologs)<-gsub("#","",names(gene_orthologs))

    # convert to character
    gene_orthologs<-gene_orthologs[,lapply(.SD,as.character),.SDcols=names(gene_orthologs)]

    # return the table
    return(gene_orthologs)
}
