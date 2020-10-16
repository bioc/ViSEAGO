#' @title Enriched GO terms intersections plot.
#' @description This method allows to visualize GO terms intersections between results of enrichment tests.
#' @importFrom data.table data.table .SD fwrite
#' @importFrom UpSetR upset
#' @family enrich_GO_terms
#' @family visualization
#' @param object an \code{\link{enrich_GO_terms-class}} or \code{\link{GO_clusters-class}} objects.
#' @param file output file name (default to "./upset.xls")
#' @details
#' This function displays and print the intersections of enriched GO terms (p<0.01) between all results provided by \code{\link{enrich_GO_terms-class}}
#' or \code{\link{GO_clusters-class}} objects. The intersections are shown in an upset plot and printed in a table.
#' @return print table and upset.
#' @include enrich_GO_terms.R GO_clusters.R
#' @examples
#' # load example object
#' data(
#'     myGOs,
#'     package="ViSEAGO"
#' )
#'
#' # print upset
#' ViSEAGO::Upset(myGOs)
#' @name Upset
#' @rdname Upset-methods
#' @exportMethod Upset
setGeneric(name="Upset",def=function(object,file="./upset.xls") {standardGeneric("Upset")})

#' @rdname Upset-methods
#' @aliases Upset
setMethod(
    "Upset",
    signature="ANY",
    definition=function(object,file){

        # check the class
        if(!is(object,"enrich_GO_terms") & !is(object,"GO_SS") & !is(object,"GO_clusters")){
            stop("object must be enrich_GO_terms, GO_SS, or GO_clusters class objects")
        }

        # extract data.table from enrich_GO_terms or GO_clusters class object.
        if(is(object,"enrich_GO_terms")){

            # extract data.table from enrich_GO_terms class object
            Data<-slot(object,"data")

        }else{

            # extract data.table from GO_clusters class object
            Data<-slot(slot(object,"enrich_GOs"),"data")
        }

        # keep only GOterms and pvalues by condition
        Data<-Data[,grep("GO\\.ID|\\.pvalue",names(Data)),with=FALSE]

        # remove pvalues in columns names
        names(Data)<-gsub("\\.pvalue","",names(Data))

        # build binary matrix for Upset graph
        Data<-data.table(
            GO.ID=Data[,"GO.ID",with=FALSE],
            Data[,lapply(.SD,function(x){val=x<0.01;x[val]<-1;x[!val]<-0;x}),.SDcols=2:ncol(Data)]
        )

        # draw upset at screen
        print(
            upset(
                Data,
                sets=rev(names(Data)[-1]),
                keep.order = TRUE,
                text.scale=2
            )
        )

        # print image to file
        png(sub("\\.[^\\.]+$",".png",file))
            print(
                upset(
                    Data,
                    sets=rev(names(Data)[-1]),
                    keep.order = TRUE,
                    text.scale=2
                )
            )
        dev.off()
        
        # build list vectors of significant GO.ID by conditions
        setlist<-lapply(2:ncol(Data),function(x){

            # extract significant GO terms
            Data$GO.ID[Data[,x,with=FALSE]==1]
        })
        names(setlist)<-names(Data)[-1]

        # use ViSEAGO internal intersetions function
        OLlist=overLapper(setlist)

        # return or print
        if(is.null(file)){

            # return
            OLlist

        }else{

            # build matrix
            OLexport <-as.matrix(
                unlist(vapply(OLlist, paste, collapse=";",""))
            )

            # convert in data.table
                OLexport<-data.table(
                combs=row.names(OLexport),
                length=lengths(OLlist),
                GOterms=OLexport[,1]
            )

            # write the file
            fwrite(
                OLexport,
                file=file,
                quote=FALSE,
                sep="\t"
            )
        }
    }
)
