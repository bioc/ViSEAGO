#' @title Display an interactive or static table.
#' @description This method is used to display or print the table for \code{\link{enrich_GO_terms-class}} or \code{\link{GO_clusters-class}} objects.
#' @importFrom methods setGeneric setMethod slot is
#' @importFrom DT datatable
#' @importFrom htmlwidgets JS
#' @family enrich_GO_terms
#' @family GO_clusters
#' @family visualization
#' @param object an \code{\link{enrich_GO_terms-class}} object from \code{\link{merge_enrich_terms}},
#' or \code{\link{GO_clusters-class}} object from \code{\link{GOterms_heatmap}}.
#' @param file table output file name (default to NULL).
#' @details This method displays an interactive table (if file=NULL) from \code{\link{enrich_GO_terms-class}} or \code{\link{GO_clusters-class}} objects.\cr
#' The table could be printed by setting file argument.
#' @references
#' Yihui Xie (2016). DT: A Wrapper of the JavaScript Library 'DataTables'. R package version 0.2. https://CRAN.R-project.org/package=DT
#' @examples
#' ###################
#' # load example object
#' utils::data(
#'  BP_sResults,
#'  package="ViSEAGO"
#' )
#'
#' ###################
#' # display merge_enrich_terms output
#' ViSEAGO::show_table(BP_sResults)
#'
#' ###################
#' # print merge_enrich_terms output
#' ViSEAGO::show_table(
#'  BP_sResults,
#'  "./data/output/BP_sResults.txt"
#' )
#'
#' \dontrun{
#' ###################
#' # display table of GO_clusters-class object
#' ViSEAGO::show_table(Wang_clusters_wardD2)
#'
#' ###################
#' # print table of GO_clusters-class object
#' ViSEAGO::show_table(
#'  Wang_clusters_wardD2,
#'  "./data/output/Wang_clusters_wardD2.txt"
#' )
#' }
#' @exportMethod show_table
setGeneric(name="show_table",def=function(object,file=NULL) {standardGeneric("show_table")})

setMethod("show_table",definition=function(object,file) {

  ###################
  # load data.table
  base::require("data.table")

  ###################
  # check the class
  if(!base::class(object)%in%c("enrich_GO_terms","GO_clusters")){
   base::stop("object must be ViSEAGO::merge_enrich_terms() or from ViSEAGO::GOterms_heatmap()")
  }

  ###################
  # remove gene identifiants for printing
  if(methods::is(object,"enrich_GO_terms")){
    data<-methods::slot(object,"data")
  }else{
    data<-methods::slot(
      methods::slot(
        object,
        "enrich_GOs"
        ),
      "data"
    )
  }

  ###################
  # if interactive mode
  if(base::is.null(file)){

    ###################
    # remove genes columns for interactive mode
    data<-data[,-grep("genes",base::names(data)),with=F]

    ###################
    # add link to amigo Gene Ontology
    #data[,`:=`(GO.ID=base::paste('<a href="http://amigo.geneontology.org/amigo/term/',GO.ID,'">',GO.ID,'</a>',sep=""))]
    data$GO.ID<-base::paste('<a href="http://amigo.geneontology.org/amigo/term/',data$GO.ID,'">',data$GO.ID,'</a>',sep="")

    ###################
    # remove . in header
    names(data)<-gsub("\\."," ",names(data))

    ###################
    # convert NA IN characters
    data[base::is.na(data)]<-"NA"

    ###################
    # create a datatable
    DT::datatable(data,

      ###################
      # escape the first column
      width=650,

      ###################
      # escape the first column
      escape=1,

      ###################
      # filter column
      filter ='top',

      ###################
      # extensions to use
      extensions =c('Buttons','FixedColumns','Scroller'),

      ###################
      # DT options to use
      options = list(

        ###################
        # DT options to use
        columnDefs = base::list(
          base::list(

            ###################
            # targets columns
            targets=if(methods::is(object,"enrich_GO_terms")){
              base::c(2:3,base::grep("genes",base::names(data)))
            }else{
              base::c(4:5,base::grep("genes",base::names(data)))
            },

            ###################
            # show the fist characters
            render =htmlwidgets::JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data.length > 50 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;",
              "}"
            )
          )
        ),

        ###################
        # the dom
        dom='Bfltipr',

        ###################
        # regex with highlight
        search = base::list(regex = TRUE, caseInsensitive = FALSE),
        searchHighlight = TRUE,

        ###################
        # x axis scroller
        scrollX = TRUE,
        fixedColumns =list(leftColumns = 1),

        ###################
        # y scroller
        deferRender = TRUE,
        scrollY = 400,
        scroller = TRUE,

        ###################
        # button download
        buttons=base::I('colvis')
      )
    )
  }else{

    ###################
    # write the table
    data.table::fwrite(data,file=file,na="NA",sep="\t")
  }
})
