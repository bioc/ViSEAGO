#' @title Store available organisms GO annotations at EntrezGene.
#' @description  Store the available species and current GO annotations from the
#' \href{https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz}{gene2go.gz} file avalable at NCBI \href{https://ftp.ncbi.nih.gov/gene/DATA}{EntrezGene ftp}.
#' @importFrom data.table data.table fread rbindlist :=
#' @family genomic_ressource
#' @details This function downloads the \href{https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz}{gene2go.gz} file from
#' \href{https://ftp.ncbi.nih.gov/gene/DATA}{EntrezGene ftp} which contains available organisms (taxid) with the corresponding GO annotations.
#' @return a  \code{\link{genomic_ressource-class}} object required by \code{\link{annotate}}.
#' @references
#' Matt Dowle and Arun Srinivasan (2017). data.table: Extension of `data.frame`. R package version 1.10.4. https://CRAN.R-project.org/package=data.table.
#'
#' Eric Sayers (2013). Entrez Programming Utilities Help.
#'
#' #' Henrik Bengtsson (2016). R.utils: Various Programming Utilities. R package version 2.5.0. https://CRAN.R-project.org/package=R.utils.
#'
#' Maglott, D, Ostell, J, Pruitt, KD, Tatusova, T (2011). Entrez Gene: gene-centered information at NCBI. Nucleic Acids Res., 39, Database issue:D52-7.
#' @include genomic_ressource.R
#' @examples
#' \dontrun{
#' ###################
#' # Download EntrezGene available organisms GO annotations
#' EntrezGene<-ViSEAGO::EntrezGene2GO()
#' }
#' @export
EntrezGene2GO=function(){

    ## download file

    # temp file
    temp<-paste(
      tempfile(),
      "gz",
      sep="."
    )

    # import Gene to Gene Ontology from NCBI Gene database
    download.file(
      "ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz",
      quiet=TRUE,
      destfile =temp
    )

    # uncompress
    gunzip(temp)

    # read the file (linux and windows)
    gene2go=fread(
      sub("\\.gz","",temp),
      verbose=FALSE,
      showProgress=FALSE
    )

    # select columns and rename
    gene2go<-unique(
      gene2go[,c(seq_len(4),8),with=FALSE]
    )
    colnames(gene2go)<-c("taxid","gene_id","GOID","evidence","category")

    # convert columns in character
    gene2go[,
      `:=`(
        taxid=as.character(gene2go$taxid),
        gene_id=as.character(gene2go$gene_id)
      )
    ]

    ## taxonomy

    # extract taxonomique informations
    taxon<-taxonomy(gene2go$taxid)

    ###################
    # return data genomic_ressource class object
    new(
        "genomic_ressource",
        db="EntrezGene",
        stamp=as.character(Sys.time()),
        data=gene2go,
        organisms=taxon
    )
}
