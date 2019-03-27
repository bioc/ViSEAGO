#' @title enrich_GO_terms class object definition.
#' @description This class is invoked by \code{\link{merge_enrich_terms}} method in order to store the merged \code{data.table} and associated metadata.
#' @importFrom methods setClass slot
#' @family enrich_GO_terms
#' @slot same_genes_background logical.
#' @slot input a list containing named elements. Each element must contain the name of \code{\link[topGO]{topGOdata-class}}
#' object created by \code{\link{create_topGOdata}} method and the associated  \code{\link[topGO]{topGOresult-class}}
#' object(s) to combinate (see examples in \code{\link{merge_enrich_terms}}).
#' @slot ont ontology used "MF", "BP", or "CC".
#' @slot topGO  a \code{list} with topGO objects summary informations.
#' @slot data a merged \code{data.table}  of enriched GO terms (p<0.01) in at least once with GO descriptions and statistical values.
setClass("enrich_GO_terms",
         slots=c(
           same_genes_background="logical",
           input="list",
           ont="character",
           topGO="list",
           data="data.table"
         )
)
#' @importFrom methods setMethod
setMethod("show", "enrich_GO_terms",function(object) {

  ###################
  # Extract table
  Data<-methods::slot(object,"data")

  ###################
  # Extract pvalues
  Data<-Data[,base::grep("\\.pvalue",base::names(Data)),with=F]

  ###################
  # count significant pvalues by condition
  Data<-Data[,
    lapply(.SD,function(x){base::sum(x<0.01,na.rm = T)}),
    .SDcols=1:base::ncol(Data)
  ]

  ###################
  # melt the table
  Data<-data.table::melt.data.table(
    Data,
    measure.vars=base::names(Data),
    variable.name = "conditions",
    value.name = "significant GO terms number"
  )

  ###################
  # remove .pvalue in conditions column
  Data[,conditions:=base::gsub("\\.pvalue","",conditions)]

  ###################
  # get topGO information
  topGO<-methods::slot(object,"topGO")

  ###################
  # format topGO information
  topGO<-base::sapply(base::names(topGO),function(x){

    ###################
    # topGO subset
    Data<-topGO[[x]]

    ###################
    # summery by element
    elems<-base::paste(base::sapply(base::names(Data),function(y){

      ###################
      # summery by element
      base::paste(
        base::paste("   ",y),
        "\n       ",
        base::paste(
          base::paste(
            base::names(Data[[y]]),
            base::sub("\n.+$","",base::unlist(Data[[y]])),
          sep=": "),
          collapse="\n        "
        )
      )
    }),collapse="\n  ")

    ###################
    # summery by element
    base::paste(base::paste(x,elems,sep ="\n  "),"\n ")
  })

  ###################
  # cat some text
  base::cat("- object class: enrich_GO_terms",
    "\n- ontology: ",methods::slot(object,"ont"),
    "\n- input:\n        ", paste(paste(base::names(object@input),
    base::sapply(methods::slot(object,"input"),function(x){base::paste(x,collapse=", ")}),sep=": "),collapse="\n        "),
    "\n- topGO summary:\n ", topGO,
    "\n- enrich GOs data.table (p<0.01 in at least one list): ",base::nrow(methods::slot(object,"data"))," GO terms of ",base::nrow(Data)," conditions.",
    base::paste("\n       ",Data$conditions,":",Data$`significant GO terms number`,"terms"),sep=""
  )
})
