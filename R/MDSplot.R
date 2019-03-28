#' @title Multi Dimensional Scale (MDS) plot
#' @description Generate a Multi Dimensional Scale (MDS) plot from distance objects.
#' @importFrom data.table data.table
#' @importFrom methods setGeneric setMethod slot is
#' @importFrom dendextend get_leaves_attr
#' @importFrom plotly plot_ly add_markers layout
#' @family GO_terms GO_clusters semantic_similarity visualization
#' @param object a \code{\link{GO_SS-class}} or \code{\link{GO_clusters-class}} objects from distances computed with \code{\link{compute_SS_distances}}.
#' @param show_clusters boolean (default to FALSE)
#' @param file static image output file name (default to NULL).
#' @details This method build and display the javascript MDSplot (if \code{file}=NULL) from \code{\link{GO_SS-class}} or \code{\link{GO_clusters-class}}
#' objects.\cr
#' A static png image could be printed by setting \code{file} argument.
#' @return a MDS plot.
#' @include GO_SS.R GO_clusters.R
#' @examples
#' ###################
#' # load objects
#' utils::data(
#'  list=base::c("myGENE2GO","BP_sResults"),
#'  package="ViSEAGO"
#' )
#'
#' ###################
#' # initialyse object for compute GO Semantic Similarity
#' myGOs<-ViSEAGO::build_GO_SS(
#'  gene2GO=myGENE2GO,
#'  enrich_GO_terms=BP_sResults
#' )
#'
#' ###################
#' # build MDS plot for a GO_SS-class distance object
#' ViSEAGO::MDSplot(myGOs)
#'
#' ##################
#' # GOtermsHeatmap with default parameters
#' Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
#'  myGOs,
#'  showIC=TRUE,
#'  showGOlabels=TRUE,
#'  GO.tree=base::list(
#'   tree=base::list(
#'    distance="Wang",
#'    aggreg.method="ward.D2",
#'    rotate=NULL
#'   ),
#'   cut=base::list(
#'    dynamic=base::list(
#'     pamStage=TRUE,
#'     pamRespectsDendro=TRUE,
#'     deepSplit=2,
#'     minClusterSize =2
#'    )
#'   )
#'  ),
#'  samples.tree=NULL
#' )
#'
#' ###################
#' # build MDS plot for a GO_clusters-class distance object, highlighting GO terms clusters.
#' ViSEAGO::MDSplot(Wang_clusters_wardD2)
#'
#' ###################
#' # compute clusters of GO terms Semantic Similarity distances
#' Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(
#'  Wang_clusters_wardD2,
#'  distance="BMA"
#' )
#'
#' ###################
#' # build MDS plot for a GO_clusters-class distance object, highlighting GO groups clusters.
#' ViSEAGO::MDSplot(
#'  Wang_clusters_wardD2,
#'  show_clusters=TRUE
#' )
#' @exportMethod MDSplot
setGeneric(name="MDSplot",def=function(object,show_clusters=FALSE,file=NULL) {standardGeneric("MDSplot")})

setMethod("MDSplot",definition=function(object,show_clusters,file) {

  #################
  # check class
  if(!base::class(object)%in%c("GO_SS","GO_clusters")){
    base::stop("object must come from ViSEAGO::compute_SS_distances() or ViSEAGO::GOterms_heatmap()")
  }

  #################
  # check class
  if(methods::is(object,"GO_SS") & show_clusters==T){
    base::stop("show_clusters is only available for GO_clusters class after clusters SS distance calculation with ViSEAGO::compute_SS_distances()")
  }

  #################
  # for SS_dist from object
  if(show_clusters==FALSE){

    #################
    # import SS_dist from object
    d=methods::slot(object,"terms_dist")

    #################
    # if empty
    if(base::length(d)==0){
      stop("Please use GO_terms class object with computed SS distances using ViSEAGO::compute_SS_distance()")
    }

    #################
    # measures
    measures=names(d)

    #################
    # MDS
    #################

      #################
      # run MDS
      res.mds=lapply(measures,function(x){

        #################
        # run MDS
        res.mds <- stats::cmdscale(d[[x]], eig =T, k = 2)

        #################
        # extract point values
        res.mds<-res.mds$points

        #################
        # convert to data.table
        data.table::data.table(
          GO.ID=base::attr(res.mds,"dimnames")[[1]],
          Dim.1=res.mds[,1],
          Dim.2=res.mds[,2],
          measure=x
        )
      })

      #################
      # bind columns
      res.mds=data.table::rbindlist(res.mds)

      #################
      # for GO_SS
      if(methods::is(object,"GO_SS")){

        #################
        # GO names
        res.mds=data.table::data.table(
          res.mds,term=methods::slot(
            methods::slot(
              object,
              "enrich_GOs"
            ),
            "data"
          )$term
        )

        #################
        # add levels to measures
        res.mds$measure<-factor(res.mds$measure,levels=unique(res.mds$measure))

      }else{

        #################
        # add clusters and GO names
        res.mds=merge(
          res.mds,
          methods::slot(
            methods::slot(
              object,
              "enrich_GOs"
            ),
            "data"
          )[,.(GO.cluster,GO.ID,term)],
          by="GO.ID",
          sort=F
          )

        #################
        # ordering by cluster
        data.table::setorder(res.mds,"GO.cluster")

        #################
        # add levels to measures
        res.mds$GO.cluster<-factor(res.mds$GO.cluster,levels=unique(res.mds$GO.cluster))
      }
  }else{

    #################
    # import SS_dist from object
    d=methods::slot(object,"clusters_dist")

    #################
    # if empty
    if(base::length(d)==0){
      stop("Please use GO_clusters class object with computed clusters distances using ViSEAGO::compute_SS_distance()")
    }

    #################
    # measures
    measures=names(d)

    #################
    # MDS
    #################

    #################
    # run MDS
    res.mds=lapply(measures,function(x){

      #################
      # run MDS
      res.mds <- stats::cmdscale(d[[x]], eig = T, k = 2)

      #################
      # extract point values
      res.mds<-res.mds$points

      #################
      # convert to data.table
      data.table::data.table(
        GO.cluster=base::attr(res.mds,"dimnames")[[1]],
        Dim.1=res.mds[,1],
        Dim.2=res.mds[,2],
        measure=x
      )
    })

    #################
    # bind columns
    res.mds=data.table::rbindlist(res.mds)

    #################
    # custom text
    res.mds[,`:=`(text=GO.cluster,GO.cluster=base::gsub("_.+$","",GO.cluster))]

    #################
    # add levels to measures
    res.mds$measure<-base::factor(res.mds$measure,levels=base::unique(res.mds$measure))

    #################
    # add levels to GO.cluster
    res.mds$GO.cluster<-base::factor(res.mds$GO.cluster,levels=base::unique(res.mds$GO.cluster))

    #################
    # add GO.ID for GO
    res.mds[,text:=base::paste("GO.cluster:",text)]

    #################
    # add GO.ID for GO
    res.mds[,text:=base::gsub("_GO","<br>GO.ID: GO",text)]

    #################
    # add GO.name for term definition
    res.mds[,text:=base::gsub("_"," <br>GO.name: ",text)]
  }

  #################
  # plot
  #################

    ################
    # mds coordinates
    values=base::unlist(res.mds$Dim.1,res.mds$Dim.2)

    ################
    # init graph
    p<-plotly::plot_ly()

    #################
    # for GO_SS
    if(methods::is(object,"GO_SS")){

      ################
      # add trace to plot by measure
      for(x in base::seq_len(length(measures))){

        ################
        # default visualization
        if(x==1){visible=TRUE}else{visible=FALSE}

        ################
        # create trace
        p<-plotly::add_markers(
          p,
          data=res.mds[measure==measures[x]],
          x=~Dim.1,
          y=~Dim.2,
         name=measures[x],
         text = ~paste('GO.ID:',GO.ID,'<br>GO.name:',term),
         showlegend=F,
          marker =base::list(
            size =20,
            opacity = 0.4,
            color="royalblue"
          ),
         visible=visible
        )
      }

      #################
      # add custom layout with dropdown menu
      p<-plotly::layout(p,

        #################
        # add title
        title="MultiDimensional Scaling plot",

        #################
        # increase font size
        font=base::list(size=14),

        #################
        # add axis legends
        xaxis=list(title="Dimension 1"),
        yaxis=list(title="Dimension 2"),

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

            buttons = base::lapply(base::seq_len(length(measures)),function(x){

              #################
              # init visibility to all FASLE
              values=base::rep(FALSE,length(measures))

              #################
              # turn to true for each measures
              values[x]<-TRUE

              #################
              # output button
              base::list(method ="restyle",
              args = base::list("visible",base::as.list(values)),
              label=measures[x])
            })
          )
        )
      )
    }

    #################
    # for GO_SS
    if(methods::is(object,"GO_clusters")){

      ###################
      # color labels of cutting tree
      colors= base::unique(
        dendextend::get_leaves_attr(methods::slot(object,"dendrograms")$GO,"edgePar")
      )

      ###################
      # extract the GO. cluster column
      count<-methods::slot(
        methods::slot(
          object,
          "enrich_GOs"
        ),
        "data"
      )[,.(GO.cluster)]

      ###################
      # count the number of terms by clusters
      count<-count[,.(nb=.N),by=GO.cluster]

      ###################
      # add count to res.mds
      res.mds=merge(res.mds,count,by="GO.cluster",sort=F)

      ###################
      # for terms
      if(show_clusters==F){

        ################
        # create trace
        p<-plotly::add_markers(
          p,
          data=res.mds,
          x=~Dim.1,
          y=~Dim.2,
          color=~GO.cluster,
          text = ~paste(
            "cluster:",
            GO.cluster,
            "<br>GO.ID:",
            GO.ID,
            "<br>GO.name:",
            term
          ),
          showlegend=T,
          colors=colors,
          marker =base::list(
            size =20,
            opacity = 0.4
          )
        )

        #################
        # add custom layout with dropdown menu
        p<-plotly::layout(p,

          #################
          # add title
          title=paste(measures,"distance MultiDimensional Scaling plot"),

          #################
          # increase font size
          font=base::list(size=14),

          #################
          # add axis legends
          xaxis=list(title="Dimension 1"),
          yaxis=list(title="Dimension 2")
        )
      }else{

        ###################
        # add count to text
        res.mds[,`:=`(
        text=paste(base::sub("<br>.+$","",text),"<br>GO.count:",
        nb,base::sub("^.+<br>GO.ID","<br>GO.ID",text)))]

        ################
        # add trace to plot by measure
        for(x in base::seq_len(base::length(measures))){

          ################
          # default visualization
          if(x==1){visible=T}else{visible=F}

          ################
          # create trace
          p<-plotly::add_markers(
            p,
            data=res.mds[measure==measures[x]],
            x=~Dim.1,
            y=~Dim.2,
            name=measures[x],
            text = ~text,
            showlegend=F,
            sizes=base::c(20,50),
            size=~nb,
            marker =base::list(
              sizemode = 'diameter',
              opacity = 0.4,
              line=base::list(color=colors),
              color=colors
            ),
            visible=visible
          )
        }

        #################
        # add custom layout with dropdown menu
        p<-plotly::layout(p,

          #################
          # add title
          title="MultiDimensional Scaling plot",

          #################
          # increase font size
          font=base::list(size=14),

          #################
          # add axis legends
          xaxis=list(title="Dimension 1"),
          yaxis=list(title="Dimension 2"),

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

              buttons = base::lapply(base::seq_len(length(measures)),function(x){

                #################
                # init visibility to all FASLE
                values=base::rep(FALSE,base::length(measures))

                #################
                # turn to true for each measures
                values[x]<-TRUE

                #################
                # output button
                base::list(method ="restyle",
                args = base::list("visible",base::as.list(values)),
                label=measures[x])
              })
            )
          )
        )
      }
    }

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
})
