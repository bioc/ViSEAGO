#' @title Barplot for the count of GO terms.
#' @description This method displays in barplot the count of GO terms splitted in two categories (significant or not)
#' for each result of GO enrichment tests.
#' @importFrom data.table data.table melt := .SD
#' @importFrom plotly plot_ly layout export
#' @importFrom methods is slot
#' @family GO_terms
#' @family visualization
#' @param object an \code{\link{enrich_GO_terms-class}} object from \code{\link{merge_enrich_terms}} method.
#' @param file the name of the output file (default to NULL for interactive screen display).
#' @details This method displays an interactive barplot, using \pkg{plotly} package, from a \code{\link{merge_enrich_terms}}
#' output object.\cr
#' A static image (in png) could be printed by setting \code{file} argument.
#' @return a barplot.
#' @references
#' Carson Sievert, Chris Parmer, Toby Hocking, Scott Chamberlain, Karthik Ram, Marianne Corvellec and Pedro Despouy (2017).
#' plotly: Create InteractiveWeb Graphics via 'plotly.js'. R package version 4.6.0. https://CRAN.R-project.org/package=plotly.
#' @include enrich_GO_terms.R
#' @examples
#' # load object
#' utils::data(
#'  myGOs,
#'  package="ViSEAGO"
#' )
#'
#' # barplot for the count of GO terms
#' ViSEAGO::GOcount( myGOs)
#' @exportMethod GOcount
#' @name GOcount
#' @rdname GOcount-methods
#' @exportMethod GOcount
setGeneric(
    name="GOcount",
    def=function(object,file=NULL){
        standardGeneric("GOcount")
    }
)

#' @rdname GOcount-methods
#' @aliases GOcount
setMethod(
    "GOcount",
    signature="ANY",
    definition=function(object,file){

        # check class
        if(!is(object,"enrich_GO_terms") & !is(object,"GO_SS") & !is(object,"GO_clusters")){
            stop("object must be enrich_GO_terms, GO_SS, or GO_clusters class objects")
        }

        # Extract data
        if(is(object,"enrich_GO_terms")){

            # data
            data<-slot(
                object,
                "data"
            )

        }else{

            # data
            data<-slot(
                slot(
                    object,
                    "enrich_GOs"
                ),
                "data"
            )
        }

        # find pvalues columns
        pvalues<-grep("\\.pvalue",names(data))

        # check more than one pvalues
        if(length(pvalues)==1){
            stop("GOcount is required for multiples comparisons")
        }

        # count significant (or not) pvalues by condition
        Data<-data[,pvalues,with=FALSE]

        # count
        Data<-data.table(
            pvalue=c("not enriched GO terms","enriched GO terms"),
            
                
                essai<-Data[,
                    lapply(.SD,function(x){
    
                        # count
                        res<-table(x<0.01)
    
                        # count
                        if(length(res)==1){
                            
                            if(names(res)=="TRUE"){c(0,res)}else{c(res,0)}
                            
                        }else{
                            as.vector(res)
                        } 
                    }),
                    .SDcols=seq_len(ncol(Data))
                ]
        )

        # remove pvalue in names
        names(Data)<-gsub("\\.pvalue","",names(Data))

        # melt the table
        Data<-melt.data.table(
            Data,
            id.vars=names(Data)[1],
            variable.name="condition",
            value.name="count"
            )

        # assign class values
        Data[,
            `:=`(
                pvalue=factor(Data$pvalue,levels = c("enriched GO terms","not enriched GO terms")),
                condition=factor(Data$condition,levels =rev(unique(Data$condition))),
                count=as.numeric(Data$count)
            )
        ]

        # barchart
        p<-plot_ly(
            data=Data,
            type = 'bar',
            y = ~condition,
            x = ~count,
            color=~pvalue,
            colors=c("#FED8B1","#87ceeb")
        )

        # barchart layout
        p<-layout(p,
            xaxis =list(
                title = 'GO terms number'
            ),
            yaxis =list(
                title =""
            ),
             margin =list(t=100,b =100,l=150),
            title="number of significant (or not) GO terms by conditions",
            barmode = 'stack'
        )

         # return or print
        if(is.null(file)){

            # return the plot
            p

        }else{

            # print heatmap
            orca(p,file=file)
        }
    }
)
