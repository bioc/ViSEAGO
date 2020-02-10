#' @title Check available organisms databases at Uniprot.
#' @description Check the \href{http://www.ebi.ac.uk/GOA}{Uniprot-GOA} available organisms.
#' @importFrom data.table data.table fread
#' @family genomic_ressource
#' @details This function downloads the current_release_numbers file (ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/current_release_numbers.txt) from
#' \href{http://www.ebi.ac.uk/GOA}{Uniprot-GOA} which contains available organisms.
#' @references
#' Matt Dowle and Arun Srinivasan (2017). data.table: Extension of `data.frame`. R package version 1.10.4. https://CRAN.R-project.org/package=data.table.
#'
#' Huntley, RP, Sawford, T, Mutowo-Meullenet, P, Shypitsyna, A, Bonilla, C, Martin, MJ, O'Donovan, C (2015).
#' The GOA database: gene Ontology annotation updates for 2015. Nucleic Acids Res., 43, Database issue:D1057-63.
#' @return a  \code{\link{genomic_ressource-class}} object required by \code{\link{annotate}}.
#' @include genomic_ressource.R
#' @examples
#' # List Uniprot-GOA available organisms
#' Uniprot<-ViSEAGO::Uniprot2GO()
#' @export
Uniprot2GO=function(){

    # temp file
    temp<-tempfile()
    
    # import Gene to Gene Ontology from NCBI Gene database
    download.file(
        "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/current_release_numbers.txt",
        quiet=TRUE,
        method="wget",
        destfile =temp
    )

    # load current release
    Data<-fread(
        temp,
        verbose=FALSE,
        showProgress=FALSE
    )

    # return data
    new(
        "genomic_ressource",
        db="Uniprot-GOA",
        stamp=as.character(Sys.time()),
        data=data.table(),
        organisms=Data
    )
}
