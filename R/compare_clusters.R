#' @title Heatmap to compare partitions
#' @description Build an interactive heatmap of the common GO terms frequency between several partitions.
#' @importFrom data.table data.table
#' @importFrom methods setGeneric setMethod slot
#' @importFrom plotly plot_ly layout
#' @family GO_clusters
#' @param clusters  a \code{list} of named \code{\link{GO_clusters-class}} objects,
#' from \code{\link{GOterms_heatmap}} or \code{\link{GOclusters_heatmap}} methods.
#' @details
#' Build an interactive heatmap of common GO terms frequency between partitions from several
#' \code{\link{GO_clusters-class}} objects.
#' @return an interactive javascript heatmap.
#' @references
#'  Carson Sievert, Chris Parmer, Toby Hocking, Scott Chamberlain, Karthik Ram, Marianne Corvellec and Pedro Despouy (2017).
#'  plotly: Create Interactive Web Graphics via 'plotly.js'. R package version 4.6.0. https://CRAN.R-project.org/package=plotly
#' @include GO_clusters.R
#' @examples
#' ###################
#' # load example object
#' utils::data(
#'  myGOs,
#'  package="ViSEAGO"
#' )
#'
#' \dontrun{
#' ###################
#' # compute Semantic Similarity (SS)
#' myGOs<-ViSEAGO::compute_SS_distances(
#'  myGOs,
#'  distance=base::c("Resnik","Lin","Rel","Jiang","Wang")
#' )
#'
#' ##################
#' # Resnik distance GO terms heatmap
#' Resnik_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#'   myGOs,
#'   showIC=TRUE,
#'   showGOlabels=TRUE,
#'   GO.tree=base::list(
#'     tree=base::list(
#'       distance="Resnik",
#'       aggreg.method="ward.D2"
#'     ),
#'     cut=base::list(
#'       dynamic=base::list(
#'         deepSplit=2,
#'         minClusterSize =2
#'       )
#'     )
#'   ),
#'   samples.tree=NULL
#' )
#'
#' ##################
#' # Lin distance GO terms heatmap
#' Lin_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#'   myGOs,
#'   showIC=TRUE,
#'   showGOlabels=TRUE,
#'   GO.tree=base::list(
#'     tree=base::list(
#'       distance="Lin",
#'       aggreg.method="ward.D2"
#'     ),
#'     cut=base::list(
#'       dynamic=base::list(
#'         deepSplit=2,
#'         minClusterSize =2
#'       )
#'     )
#'   ),
#'   samples.tree=NULL
#' )
#'
#' ##################
#' # Rel distance GO terms heatmap
#' Rel_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#'   myGOs,
#'   showIC=TRUE,
#'   showGOlabels=TRUE,
#'   GO.tree=base::list(
#'     tree=base::list(
#'       distance="Rel",
#'       aggreg.method="ward.D2"
#'     ),
#'     cut=base::list(
#'       dynamic=base::list(
#'         deepSplit=2,
#'         minClusterSize =2
#'       )
#'     )
#'   ),
#'   samples.tree=NULL
#' )
#'
#' ##################
#' # Jiang distance GO terms heatmap
#' Jiang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#'   myGOs,
#'   showIC=TRUE,
#'   showGOlabels=TRUE,
#'   GO.tree=base::list(
#'     tree=base::list(
#'       distance="Jiang",
#'       aggreg.method="ward.D2"
#'     ),
#'     cut=base::list(
#'       dynamic=base::list(
#'         deepSplit=2,
#'         minClusterSize =2
#'       )
#'     )
#'   ),
#'   samples.tree=NULL
#' )
#'
#' ##################
#' # Wang distance GO terms heatmap
#' Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#'   myGOs,
#'   showIC=TRUE,
#'   showGOlabels=TRUE,
#'   GO.tree=base::list(
#'     tree=base::list(
#'       distance="Wang",
#'       aggreg.method="ward.D2"
#'     ),
#'     cut=base::list(
#'       dynamic=base::list(
#'         deepSplit=2,
#'         minClusterSize =2
#'       )
#'     )
#'   ),
#'   samples.tree=NULL
#' )
#' }
#'
#' ###################
#' # clusters to compare
#' clusters<-base::list(
#'  Resnik="Resnik_clusters_wardD2",
#'  Lin="Lin_clusters_wardD2",
#'  Rel="Rel_clusters_wardD2",
#'  Jiang="Jiang_clusters_wardD2",
#'  Wang="Wang_clusters_wardD2"
#' )
#'
#' \dontrun{
#' ###################
#' # clusters content comparisons
#' clusters_comp<-ViSEAGO::compare_clusters(clusters)
#' }
#' @name compare_clusters
#' @rdname compare_clusters-methods
#' @exportMethod compare_clusters
setGeneric(name="compare_clusters",def=function(clusters) {standardGeneric("compare_clusters")})


#' @rdname compare_clusters-methods
#' @aliases compare_clusters
setMethod("compare_clusters",signature ="list",definition=function(clusters) {

  ###################
  # convert GO.ID GO.clusters column in list
  ###################

    ###################
    # extract clusters and GO.ID
    Clusters<-base::unlist(recursive=F,base::lapply(clusters,function(x){

      ###################
      # load the target object
      object=base::get(x[[1]])

      #################
      # check class
      if(!methods::is(object,"GO_clusters")) base::stop("object must come from ViSEAGO::GOterms_heatmap()")

      ###################
      # extract GO.ID and correspond cluster assignation to a list
      object=utils::unstack(
        methods::slot(
          methods::slot(
            object,
            "enrich_GOs"
          ),
          "data"
        )[,base::c("GO.ID","GO.cluster"),with=FALSE]
      )

      ###################
      # names of list elements
      list.names<-base::as.numeric(base::names(object))

      ###################
      # add 0 if needed
      list.names<-sprintf(paste("%0",base::max(base::nchar(list.names)),"d",sep=""), list.names)

      ###################
      # rename list elements
      base::names(object)<-paste(base::names(x),list.names,sep="cl")

      ###################
      # return the object
      object
    }))

  ###################
  # count common terms
  ###################

    ###################
    # return a matrix cross-product from all common GO.ID count by Cluster-objects
    clusters.length<-base::lengths(Clusters)

    ###################
    # return a matrix cross-product from all common GO.ID count by Cluster-objects
    Data<-base::crossprod(base::table(utils::stack(Clusters)))

    ###################
    # assign 0 to upper part of the matrix
    Data[base::upper.tri(Data)] <-0

    ###################
    # reverse matrix row
    Data<-Data[base::nrow(Data):1,]

  ###################
  # rectangle shape definition
  ###################

    ###################
    # methods
    methods=base::unique(gsub("\\..+$","",base::colnames(Data)))

    ###################
    # rectangle shape definition
    coord<-base::unlist(recursive =F,lapply(base::seq_len(base::length(methods)),function(x){

      ###################
      # x coordinates
      col=base::grep(base::paste(methods[x],"\\.",sep=""),base::colnames(Data))

      ###################
      # y coordinates
      base::lapply(x:base::length(methods),function(y){

        ###################
        # y coordinates
        row=base::grep(base::paste(methods[y],".",sep=""),base::row.names(Data))

        ###################
        # y coordinates
        val=c(col[1]-1.5,utils::tail(col,n=1)-0.5,row[1]-1.5,utils::tail(row,n=1)-0.5)

        ###################
        # val names
        base::names(val)<-c("x0","x1","y0","y1")

        ###################
        # return val
        val
      })
    }))

    ###################
    # initiate a rect shape object
    rect <- base::list(
      type = "rect",
      fillcolor = "",
      line = base::list(color = "black",width=1,dash="dash"),
      opacity = 0.3,
      xref = "x",
      yref = "y"
    )

    ###################
    # all rectangles
    rects <-base::lapply(coord,function(x){

      ###################
      # rect dimensions
      rect[["x0"]] <- x["x0"]
      rect[["x1"]] <- x["x1"]
      rect[["y0"]] <- x["y0"]
      rect[["y1"]] <- x["y1"]

      ###################
      # rect rect
      rect
    })

  ###################
  # build plot matrix
  ###################

    ###################
    # count by line
    rowsum=base::apply(Data,1,sum)

    ###################
    # count by column
    colsum=base::apply(Data,2,sum)

    ###################
    # keep Data without all 0
    #Data<-Data[rowsum>0,colsum>0]

    ###################
    # keep Data without all 0
    count<-matrix(paste('<br>GO.common: ',Data),ncol=base::ncol(Data),byrow=F)

    ###################
    # row %
    row_prc<-base::round((Data/clusters.length[base::row.names(Data)])*100,digits=1)
    # row_prc[base::is.na(row_prc)]<-0

    ###################
    # col % hover text
    row_prc_hover.text=matrix(base::paste('<br>y%:',row_prc),ncol=base::ncol(Data))

    ###################
    # replace 0 by NA
    row_prc[row_prc==0]<-NA

    ###################
    # col %
    col_prc<-base::t(base::round((base::t(Data)/clusters.length[base::colnames(Data)])*100,digits=1))

    ###################
    # col % hover text
    col_prc_hover.text=matrix(base::paste('<br>x%:',col_prc),ncol=base::ncol(Data))

    ###################
    # replace 0 by NA
    col_prc[col_prc==0]<-NA

    ###################
    # column GOcount
    x=matrix(paste(
        'x:',
        base::rep(base::colnames(col_prc),each=base::ncol(col_prc)),
        '<br>x GO sum:',
        base::rep(clusters.length[base::colnames(col_prc)],each=base::ncol(col_prc))),
      ncol=base::ncol(col_prc)
    )

    ###################
    # row GOcount
    y=matrix(paste(
      '<br>y:',
      base::rep(base::row.names(col_prc),each=base::ncol(col_prc)),
      '<br>y GO sum:',
      base::rep(clusters.length[base::row.names(col_prc)],
      each=base::ncol(col_prc))),ncol=base::ncol(col_prc),
    byrow=TRUE)

    ###################
    # create hover text matrix
    hover.text<-base::matrix(paste(x,y,count,col_prc_hover.text,
    row_prc_hover.text),ncol=base::ncol(count),byrow=F)

    ###################
    # remove hover text with no GO in common
    hover.text[  count=="<br>GO.common:  0"]<-""

  ###################
  # interactive heatmap
  ###################

    ###################
    # initialyse p
    p<-plotly::plot_ly()

    ###################
    # Draw heatmap with y%
    p<-plotly::add_heatmap(p,x=colnames(col_prc),y=row.names(col_prc),z =col_prc,
    colorscale = "Greys",name="x identity%",hoverinfo='text',text=hover.text,visible=TRUE)

    ###################
    # Draw heatmap with y%
    p<-plotly::add_heatmap(p,x=colnames(row_prc),y=row.names(row_prc),z =row_prc,
    colorscale = "Greys",name="y identity%",hoverinfo='text',text=hover.text,visible=FALSE)

    #################
    # add custom layout with dropdown menu
    plotly::layout(p,

      #################
      # add title
      title="clusters comparisons",

      #################
      # add axis legends
      xaxis=list(title="x",tickfont=base::list(size=10)),
      yaxis=list(title="y",tickfont=base::list(size=10)),

      ##################
      # custom margin
      margin=base::list(t=30,b=100,l=100,r=0),

      #################
      # add scrolling menu for availables measures
      shapes = rects,

      #################
      # add scrolling menu for availables measures
      updatemenus = base::list(

      #################
      # measure dropdown menu
      base::list(

        #################
        # button position
        x = 0.1,
        y = 1.1,

        buttons = base::list(

          #################
          # x output button
          base::list(method ="restyle",args = base::list("visible",base::list(TRUE,FALSE)),label="x identity%"),

          #################
          # y output button
          base::list(method ="restyle",args = base::list("visible",base::list(FALSE,TRUE)),label="y indentity%")

        )
      )
    )
  )
})
