#' @title Store organisms GO annotations from custom database file.
#' @description  Store the available species and current GO annotations from a custom table file
#' @importFrom data.table data.table fread rbindlist :=
#' @family genomic_ressource
#' @param file custom GO annotation file
#' @details This function load a custom GO annotation database table with columns:
#' \describe{
#'      \item{taxid}{custom taxonomic identifiants}
#'      \item{gene_id}{custom gene identifiants}
#'      \item{GOID}{Known GO identifiants (see \code{select(GO.db,columns=columns(GO.db),keys=keys(GO.db))}}
#'      \item{evidence}{Known GO \href{http://geneontology.org/page/guide-go-evidence-codes}{evidence codes}}
#' }
#' @return a  \code{\link{genomic_ressource-class}} object required by \code{\link{annotate}}.
#' @references
#' Matt Dowle and Arun Srinivasan (2017). data.table: Extension of `data.frame`. R package version 1.10.4. https://CRAN.R-project.org/package=data.table.
#' @include genomic_ressource.R
#' @examples
#' ###################
#' # Download custom GO annotations
#' Custom<-ViSEAGO::Custom2GO(
#'     system.file(
#'         "extdata/customfile.txt",
#'         package = "ViSEAGO"
#'     )
#' )
#' @export
Custom2GO=function(file){

    # read the file (linux and windows)
    gene2go<-fread(file)

    # check columns name
    if(!all(names(gene2go)%in%c("taxid","gene_id","GOID","evidence"))){
        stop('custom annotation file columns names required: "taxid","gene_id","GOID","evidence"')
    }

    # convert columns in character
    gene2go[,
        `:=`(
            taxid=as.character(gene2go$taxid),
            gene_id=as.character(gene2go$gene_id)
        )
    ]

    # check GOID validity
    if(!all(gene2go$GOID%in%keys(GO.db))){
        stop('GOID must be a valid identifiant. see GO.db,columns=columns(GO.db),keys=keys(GO.db)')
    }

    # extract Go category
    GO<-select(GO.db,columns=c("GOID","ONTOLOGY"),keys=unique(gene2go$GOID))

    # renanme ontology column
    names(GO)[2]<-"category"

    # merge with Gene2GO
    gene2go<-merge(gene2go,GO,by="GOID",all.x=TRUE,sort=FALSE)

    ###################
    # return data genomic_ressource class object
    new(
        "genomic_ressource",
        db="Custom",
        stamp=as.character(Sys.time()),
        data=gene2go,
        organisms=data.table(taxid=unique(gene2go$taxid))
    )
}
