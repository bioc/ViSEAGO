#' @title GO_SS class object definition.
#' @description This class is invoked by \code{\link{build_GO_SS}} method in order to store \code{\link{enrich_GO_terms-class}} object, Information Content (IC),
#' and GO terms or groups distances objects based on semantic similarity.
#' @importFrom methods setClass
#' @importFrom data.table melt.data.table .SD
#' @family GO_semantic_similarity
#' @slot db should be "Bioconductor", "EntrezGene", "Ensembl", or "Uniprot" ressource name.
#' @slot stamp date of stamp or annotation version if available in \code{character}.
#' @slot organism target species GO anotation in \code{character}.
#' @slot ont used ontology with "MF", "BP", or "CC".
#' @slot topGO  \code{list} with topGO objects summary informations.
#' @slot IC Information Content (IC)
#' @slot enrich_GOs \code{\link{merge_enrich_terms}} output object (\code{\link{enrich_GO_terms-class}} object).
#' @slot terms_dist \code{list} of GO terms or groups distances objects based on semantic similarity.
#' @include enrich_GO_terms.R
setClass("GO_SS",
         slots=c(
           db="character",
           stamp ="character",
           organism="character",
           ont="character",
           topGO="list",
           IC ="numeric",
           enrich_GOs="enrich_GO_terms",
           terms_dist="list"
        )
)

#' @importFrom methods setMethod slot is
setMethod("show",signature="GO_SS",function(object){

  ###################
  # Extract table
  Data<-methods::slot(
    methods::slot(
      object,
      "enrich_GOs"
      ),
    "data"
  )

  ###################
  # Extract pvalues
  Data<-Data[,base::grep("\\.pvalue",base::names(Data)),with=F]

  ###################
  # count significant pvalues by condition
  Data<-Data[,lapply(.SD,function(x){base::sum(x<0.01)}),
    .SDcols=base::seq_len(base::ncol(Data))
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
  Data[,"conditions":=base::gsub("\\.pvalue","",Data$conditions)]

  ###################
  # get topGO information
  topGO<-methods::slot(object,"topGO")

  ###################
  # format topGO information
  topGO<-base::vapply(base::names(topGO),function(x){

    ###################
    # topGO subset
    Data<-topGO[[x]]

    ###################
    # summery by element
    elems<-base::paste(base::vapply(base::names(Data),function(y){

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
    },""),collapse="\n  ")

    ###################
    # summery by element
    base::paste(base::paste(x,elems,sep ="\n  "),"\n ")
  },"")

  ###################
  # cat some text
  base::cat("- object class: GO_SS",
    "\n- database: ",
    methods::slot(object,"db"),
    "\n- stamp/version: ",
    methods::slot(object,"stamp"),
    "\n- organism id: ",
    methods::slot(object,"organism"),
    "\n- ontology: ",
    methods::slot(object,"ont"),
    "\n- input:\n        ",
    paste(
      paste(
        base::names(
          methods::slot(
            methods::slot(
              object,
              "enrich_GOs"
              ),
            "input"
          )
        ),
      base::vapply(
        methods::slot(
          methods::slot(
            object,
            "enrich_GOs"
          ),
          "input"
        ),
        function(x){base::paste(x,collapse=", ")}
      ,""),
      sep=": "
      ),
      collapse="\n        "
    ),
    "\n- topGO summary:\n ",
    topGO,
    "\n- enrich GOs data.table: ",
    base::nrow(
      methods::slot(
        methods::slot(
          object,
          "enrich_GOs"
        ),
        "data"
      )
    ),
    " GO terms of ",
    base::nrow(Data),
    " conditions.",
    base::paste(
      "\n       ",
      Data$conditions,
      ":",
      Data$`significant GO terms number`,
      "terms"
    ),
    if(base::length(object@terms_dist)>0){
      base::paste(
        "\n- terms distances: ",
        base::paste(
          base::names(
            methods::slot(object,"terms_dist")
          ),
          collapse=", "
        )
      )
    },
    sep=""
  )
})
