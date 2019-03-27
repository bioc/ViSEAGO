#' @title Barplot for the count of GO terms.
#' @description This method displays in barplot the count of GO terms splitted in two categories (significant or not)
#' for each result of GO enrichment tests.
#' @importFrom data.table data.table melt :=
#' @importFrom methods setGeneric setMethod
#' @importFrom plotly plot_ly layout add_trace
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
#' ###################
#' # load object
#' utils::data(
#'  BP_sResults,
#'  package="ViSEAGO"
#' )
#'
#' ###################
#' # barplot for the count of GO terms
#' ViSEAGO::GOcount(BP_sResults)
#' @exportMethod GOcount
setGeneric(name="GOcount",def=function(object,file=NULL){standardGeneric("GOcount")})

setMethod("GOcount",signature="enrich_GO_terms",definition=function(object,file) {

  ###################
  # data
  data=methods::slot(object,"data")

  ###################
  # find pvalues columns
  pvalues=base::grep("\\.pvalue",base::names(data))

  ###################
  # check more than one pvalues
  if(base::length(pvalues)==1){

    base::cat("GOcount is required for multiples comparisons")

  }else{

    ###################
    # count significant (or not) pvalues by condition
    Data<-data[,pvalues,with=F]

    ###################
    # count
    Data<-data.table::data.table(
      pvalue=c("not enriched GO terms","enriched GO terms"),
      Data[,base::lapply(.SD,function(x){

        ###################
        # count
        res<-table(x<0.01)

        ###################
        #
        if(length(res)==1){

          if(base::names(res)=="TRUE"){base::c(0,res)}else{base::c(res,0)}

        }else{

          res
        }
        }),.SDcols=base::seq_len(base::ncol(Data))]
    )

    ###################
    # remove pvalue in names
    base::names(Data)<-base::gsub("\\.pvalue","",base::names(Data))

    ###################
    # melt the table
    Data<-data.table::melt.data.table(
      Data,
      id.vars=base::names(Data)[1],
      variable.name="condition",
      value.name="count"
    )

    ###################
    # assign class values
    Data[,`:=`(
      pvalue=base::factor(pvalue,levels = c("enriched GO terms","not enriched GO terms")),
      condition=base::factor(condition,levels =base::rev(base::unique(condition))),
      count=base::as.numeric(count)
    )]

    ###################
    # barchart
    p<-plotly::plot_ly(
      data=Data,
      type = 'bar',
      y = ~condition,
      x = ~count,
      color=~pvalue,
      colors=c("#FED8B1","#87ceeb")
   #   opacity=0.8
    )

    p<-plotly::layout(p,
      xaxis =base::list(
        title = 'GO terms number'
        ),
      yaxis =base::list(
        title =""
      ),
      margin =base::list(b =100,l=150),
      title="number of significant (or not) GO terms by conditions",
      barmode = 'stack'
    )

    #################
    # return or print
    if(base::is.null(file)){

      #################
      # return the plot
      p

    }else{

      ##################
      # print heatmap
      plotly::export(p,file=file)
    }
  }
})
