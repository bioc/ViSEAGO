#' @title Display an interactive or static table.
#' @description This method is used to display or print the table for \code{\link{enrich_GO_terms-class}} or \code{\link{GO_clusters-class}} objects.
#' @importFrom DT datatable
#' @importFrom data.table fwrite
#' @importFrom htmlwidgets JS
#' @family enrich_GO_terms
#' @family GO_clusters
#' @family visualization
#' @param object an \code{\link{enrich_GO_terms-class}} object from \code{\link{merge_enrich_terms}},
#' or \code{\link{GO_clusters-class}} object from \code{\link{GOterms_heatmap}}.
#' @param file table output file name (default to NULL).
#' @details This method displays an interactive table (if file=NULL) from \code{\link{enrich_GO_terms-class}} or \code{\link{GO_clusters-class}} objects.\cr
#' The table could be printed by setting file argument.
#' @return display or print table
#' @references
#' Yihui Xie (2016). DT: A Wrapper of the JavaScript Library 'DataTables'. R package version 0.2. https://CRAN.R-project.org/package=DT
#' @examples
#' ###################
#' # load example object
#' data(
#'     myGOs,
#'     package="ViSEAGO"
#' )
#'
#' ###################
#' # display merge_enrich_terms output
#' ViSEAGO::show_table(myGOs)
#'
#' ###################
#' # print merge_enrich_terms output
#' ViSEAGO::show_table(
#'     myGOs,
#'     "myGOs.txt"
#' )
#'
#' \dontrun{
#' ###################
#' # compute GO terms Semantic Similarity distances
#' myGOs<-ViSEAGO::compute_SS_distances(
#'     distance="Wang"
#' )
#'
#' ##################
#' # GOtermsHeatmap with default parameters
#' Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#'     myGOs,
#'     showIC=TRUE,
#'     showGOlabels=TRUE,
#'     GO.tree=list(
#'         tree=list(
#'             distance="Wang",
#'             aggreg.method="ward.D2",
#'             rotate=NULL
#'         ),
#'         cut=list(
#'             dynamic=list(
#'                 pamStage=TRUE,
#'                 pamRespectsDendro=TRUE,
#'                 deepSplit=2,
#'                 minClusterSize =2
#'             )
#'         )
#'     ),
#'     samples.tree=NULL
#' )
#' ###################
#' # display table of GO_clusters-class object
#' ViSEAGO::show_table(Wang_clusters_wardD2)
#'
#' ###################
#' # print table of GO_clusters-class object
#' ViSEAGO::show_table(
#'     Wang_clusters_wardD2,
#'     "Wang_clusters_wardD2.txt"
#' )
#' }
#' @name show_table
#' @rdname show_table-methods
#' @exportMethod show_table
setGeneric(
    name="show_table",
    def=function(object,file=NULL) {
        standardGeneric("show_table")
    }
)

#' @rdname show_table-methods
#' @aliases show_table
setMethod(
    "show_table",
    signature="ANY",
    definition=function(object,file){

        # check class
        if(!is(object,"enrich_GO_terms") & !is(object,"GO_SS") & !is(object,"GO_clusters")){
            stop("object must be enrich_GO_terms, GO_SS, or GO_clusters class objects")
        }

        # Extract data
        if(is(object,"enrich_GO_terms")){

            # data
            data<-slot(object,"data")

        }else{

            # data
            data<-slot(slot(object,"enrich_GOs"),"data")
        }

        # if interactive mode
        if(is.null(file)){

            # remove genes columns for interactive mode
             data<-data[,-grep("genes",names(data)),with=FALSE]

            # add link to amigo Gene Ontology
            data$GO.ID<-paste('<a href="http://amigo.geneontology.org/amigo/term/',data$GO.ID,'">',data$GO.ID,'</a>',sep="")

            # remove . in header
            names(data)<-gsub("\\."," ",names(data))

            # convert NA in characters
            data[is.na(data)]<-"NA"

            # create a datatable
            datatable(data,

                # escape the first column
                width=650,

                # escape the first column
                escape=1,

                # filter column
                filter ='top',

                # extensions to use
                extensions =c('Buttons','FixedColumns','Scroller'),

                # DT options to use
                options = list(

                    # DT options to use
                    columnDefs = list(
                        list(

                            # targets columns
                            targets=if(is(object,"enrich_GO_terms")){
                                c(2:3,grep("genes",names(data)))
                            }else{
                                c(4:5,grep("genes",names(data)))
                            },

                            # show the fist characters
                            render=JS(
                                "function(data, type, row, meta) {",
                                "return type === 'display' && data.length > 50 ?",
                                "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;",
                                "}"
                            )
                        )
                    ),

                    # the dom
                    dom='Bfltipr',

                    # regex with highlight
                    search = list(regex=TRUE,caseInsensitive=FALSE),
                    searchHighlight = TRUE,

                    # x axis scroller
                    scrollX = TRUE,
                    fixedColumns =list(leftColumns = 1),

                    # y scroller
                    deferRender = TRUE,
                    scrollY = 400,
                    scroller = TRUE,

                    # button download
                    buttons=I('colvis')
                )
            )
        }else{

            # write the table
            fwrite(data,file=file,na="NA",sep="\t")
        }
    }
)
