#' @title Enriched GO terms intersections plot.
#' @description This method allows to visualize GO terms intersections between results of enrichment tests.
#' @importFrom data.table data.table .SD
#' @importFrom methods setGeneric setMethod slot is
#' @importFrom UpSetR upset
#' @family enrich_GO_terms
#' @family visualization
#' @param object an \code{\link{enrich_GO_terms-class}} object from \code{\link{merge_enrich_terms}}
#' @param file output file name (default to "./upset.xls")
#' @details
#' This function displays the intersections of enriched GO terms (p<0.01) between all results provided by \code{\link{enrich_GO_terms-class}}
#' or \code{\link{GO_clusters-class}} objects. The intersections are shown in an upset plot and printed in a table.
#' @include enrich_GO_terms.R GO_clusters.R
#' @examples
#' ##################
#' # load object
#' utils::data(
#'  list=base::c("BP_sResults","Wang_clusters_wardD2"),
#'  package="ViSEAGO"
#' )
#'
#' ##################
#' # print upset
#' ViSEAGO::Upset(BP_sResults)
#'
#' ##################
#' # print upset
#' ViSEAGO::Upset(Wang_clusters_wardD2)
#' @exportMethod Upset
setGeneric(name="Upset",def=function(object,file="./upset.xls") {standardGeneric("Upset")})

setMethod("Upset",definition=function(object,file){

  ###################
  # check the class
  if(!base::class(object)%in%c("enrich_GO_terms","GO_clusters"))base::stop("object must be ViSEAGO::merge_enrich_terms() or from ViSEAGO::GOterms_heatmap()")

  ##################
  # extract data.table from enrich_GO_terms or GO_clusters class object.
  if(methods::is(object,"enrich_GO_terms")){

    ##################
    # extract data.table from enrich_GO_terms class object
    Data<-methods::slot(object,"data")

  }else{

    ##################
    # extract data.table from GO_clusters class object
    Data<-methods::slot(
      methods::slot(
        object,
        "enrich_GOs"
        ),
      "data"
    )
  }

  ##################
  # keep only GOterms and pvalues by condition
  Data<-Data[,base::grep("GO\\.ID|\\.pvalue",base::names(Data)),with=F]

  ##################
  # remove pvalues in columns names
  base::names(Data)<-base::gsub("\\.pvalue","",base::names(Data))

  ##################
  # build binary matrix for Upset graph
  Data<-data.table::data.table(
    GO.ID=Data[,GO.ID],
    Data[,base::lapply(.SD,function(x){
      val=x<0.01;x[val]<-1;x[!val]<-0;x}),
      .SDcols=2:base::ncol(Data)]
  )

  ##################
  # draw upset
  UpSetR::upset(Data,sets=base::rev(names(Data)[-1]),keep.order = T,text.scale=3)

  ##################
  # build list vectors of significant GO.ID by conditions
  setlist<-base::lapply(2:base::ncol(Data),function(x){

    ##################
    # extract significant GO terms
    Data$GO.ID[Data[,x,with=F]==1]
  })
  base::names(setlist)<-base::names(Data)[-1]

  ##################
  # use ViSEAGO internal intersetions function
  OLlist=ViSEAGO::overLapper(setlist)

  ##################
  # return or print
  if(base::is.null(file)){

    ##################
    # return
    OLlist

  }else{

    ##################
    # build matrix
    OLexport <-base::as.matrix(
      base::unlist(base::vapply(OLlist, base::paste, collapse=";",""))
    )

    ##################
    # convert in data.table
    OLexport<-data.table::data.table(
      combs=base::row.names(OLexport),
      length=lengths(OLlist),
      GOterms=OLexport[,1]
    )

    ##################
    # write the file
    utils::write.table(OLexport, file=file, row.names=F, quote=F, sep="\t")
  }
})
