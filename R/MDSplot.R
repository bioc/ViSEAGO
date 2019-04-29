#' @title Multi Dimensional Scale (MDS) plot
#' @description Generate a Multi Dimensional Scale (MDS) plot from distance objects.
#' @importFrom data.table data.table rbindlist setorder := .N
#' @importFrom dendextend get_leaves_attr
#' @importFrom plotly plot_ly add_markers layout export
#' @family GO_terms GO_clusters semantic_similarity visualization
#' @param object a \code{\link{GO_SS-class}} or \code{\link{GO_clusters-class}} objects from distances computed with \code{\link{compute_SS_distances}}.
#' @param type could be "GOterms" to display GOterms MDSplot, or "GOclusters" to display GOclusters MDSplot.
#' @param file static image output file name (default to NULL).
#' @details This method build and display the javascript MDSplot (if \code{file}=NULL) from \code{\link{GO_SS-class}} or \code{\link{GO_clusters-class}}
#' objects.\cr
#' A static png image could be printed by setting \code{file} argument.
#' @return a MDS plot.
#' @include GO_SS.R GO_clusters.R
#' @examples
#' ###################
#' # load data example
#' utils::data(
#'  myGOs,
#'  package="ViSEAGO"
#' )
#' \dontrun{
#' ###################
#' # compute GO terms Semantic Similarity distances
#' myGOs<-ViSEAGO::compute_SS_distances(
#'     myGOs,
#'     distance="Wang"
#' )
#'
#' ###################
#' # build MDS plot for a GO_SS-class distance object
#' ViSEAGO::MDSplot(myGOs)
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
#'
#' ###################
#' # build MDS plot for a GO_clusters-class distance object, highlighting GO terms clusters.
#' ViSEAGO::MDSplot(
#'     Wang_clusters_wardD2,
#'     "GOterms"
#' )
#'
#' ###################
#' # compute clusters of GO terms Semantic Similarity distances
#' Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(
#'     Wang_clusters_wardD2,
#'     distance="BMA"
#' )
#'
#' ###################
#' # GOclusters heatmap
#' Wang_clusters_wardD2<-ViSEAGO::GOclusters_heatmap(
#'     Wang_clusters_wardD2,
#'     tree=list(
#'         distance="BMA",
#'         aggreg.method="ward.D2",
#'         rotate=NULL
#'     )
#' )
#'
#' ###################
#' # build MDS plot for a GO_clusters-class distance object, highlighting GO groups clusters.
#' ViSEAGO::MDSplot(
#'     Wang_clusters_wardD2,
#'     "GOclusters"
#' )
#' }
#' @exportMethod MDSplot
#' @name MDSplot
#' @rdname MDSplot-methods
#' @exportMethod MDSplot
setGeneric(
    name="MDSplot",
    def=function(object,type="GOterms",file=NULL){
        standardGeneric("MDSplot")
    }
)

#' @rdname MDSplot-methods
#' @aliases MDSplot
setMethod(
    "MDSplot",
    signature="ANY",
    definition=function(object,type,file){

        # check class
        if(!is(object,"GO_SS") & !is(object,"GO_clusters")){
            stop("object must come from ViSEAGO::compute_SS_distances() or ViSEAGO::GOterms_heatmap()")
        }

        # check class
        if(is(object,"GO_SS") & type=="GOclusters"){
            stop("show_clusters is only available for GO_clusters class after clusters SS distance calculation with ViSEAGO::compute_SS_distances()")
        }

        # check type argument
        type=match.arg(type,c("GOterms","GOclusters"))

        # for SS_dist from object
        if(type=="GOterms"){

            # import SS_dist from object
            d=slot(object,"terms_dist")

            # if empty
            if(length(d)==0){
                stop("Please use GO_terms class object with computed SS distances using ViSEAGO::compute_SS_distance()")
            }

            # measures
            measures=names(d)

            ## MDS

            # run MDS
            res.mds=lapply(measures,function(x){

                # run MDS
                res.mds <-cmdscale(d[[x]], eig =T, k = 2)

                # extract point values
                res.mds<-res.mds$points

                # convert to data.table
                data.table(
                    GO.ID=attr(res.mds,"dimnames")[[1]],
                    Dim.1=res.mds[,1],
                    Dim.2=res.mds[,2],
                    measure=x
                )
            })

            # bind columns
            res.mds=rbindlist(res.mds)

            # for GO_SS
            if(is(object,"GO_SS")){

                # GO names
                res.mds=data.table(
                    res.mds,
                    term=slot(slot(object,"enrich_GOs"),"data")$term
                )

                # add levels to measures
                res.mds$measure<-factor(
                    res.mds$measure,
                    levels=unique(res.mds$measure)
                )

            }else{

                # add clusters and GO names
                res.mds=merge(
                    res.mds,
                    slot(slot(object,"enrich_GOs"),"data")[,c("GO.cluster","GO.ID","term"),with=FALSE],
                    by="GO.ID",
                    sort=FALSE
                )

                # ordering by cluster
                setorder(res.mds,"GO.cluster")

                # add levels to measures
                res.mds$GO.cluster<-factor(
                    res.mds$GO.cluster,
                    levels=unique(res.mds$GO.cluster)
                )
            }
        }else{

            # import SS_dist from object
            d=slot(object,"clusters_dist")

            # if empty
            if(length(d)==0){
                stop("Please use GO_clusters class object with computed clusters distances using ViSEAGO::compute_SS_distance()")
            }

            # measures
            measures=names(d)

            ## MDS

            # run MDS
            res.mds=lapply(measures,function(x){

                # run MDS
                res.mds <-cmdscale(d[[x]], eig = T, k = 2)

                # extract point values
                res.mds<-res.mds$points

                # convert to data.table
                data.table(
                    GO.cluster=attr(res.mds,"dimnames")[[1]],
                    Dim.1=res.mds[,1],
                    Dim.2=res.mds[,2],
                    measure=x
                )
            })

            # bind columns
            res.mds=rbindlist(res.mds)

            # custom text
            res.mds[,`:=`(text=res.mds$GO.cluster,GO.cluster=gsub("_.+$","",res.mds$GO.cluster))]

            # add levels to measures
            res.mds$measure<-factor(
                res.mds$measure,
                levels=unique(res.mds$measure)
            )

            # add levels to GO.cluster
            res.mds$GO.cluster<-factor(
                res.mds$GO.cluster,
                levels=unique(res.mds$GO.cluster)
            )

            # add GO.ID for GO
            res.mds[,text:=paste("GO.cluster:",text)]

            # add GO.ID for GO
            res.mds[,text:=gsub("_GO","<br>GO.ID: GO",text)]

            # add GO.name for term definition
            res.mds[,text:=gsub("_"," <br>GO.name: ",text)]
        }

        ## plot

        # mds coordinates
        values=unlist(res.mds$Dim.1,res.mds$Dim.2)

        # init graph
        p<-plot_ly()

        # for GO_SS
        if(is(object,"GO_SS")){

            # add trace to plot by measure
            for(x in seq_along(measures)){

                # default visualization
                if(x==1){visible=TRUE}else{visible=FALSE}

                # create trace
                p<-add_markers(
                    p,
                    data=res.mds[measures[x],on="measure"],
                    x=~Dim.1,
                    y=~Dim.2,
                    name=measures[x],
                    text=~paste('GO.ID:',GO.ID,'<br>GO.name:',term),
                    showlegend=FALSE,
                    marker =list(
                        size =20,
                        opacity = 0.4,
                        color="royalblue"
                    ),
                    visible=visible
                )
            }

            # add custom layout with dropdown menu
            p<-layout(
                p,

                # add title
                title="MultiDimensional Scaling plot",

                # increase font size
                font=list(size=14),

                # add axis legends
                xaxis=list(title="Dimension 1"),
                yaxis=list(title="Dimension 2"),

                # add scrolling menu for availables measures
                updatemenus = list(

                    # measure dropdown menu
                    list(

                        # button position
                        x = 0.1,
                        y = 1.1,

                        buttons = lapply(seq_along(measures),function(x){

                            # init visibility to all FASLE
                            values=rep(FALSE,length(measures))

                            # turn to true for each measures
                            values[x]<-TRUE

                            # output button
                            list(
                                method="restyle",
                                args=list("visible",as.list(values)),
                                label=measures[x]
                            )
                        })
                    )
                )
            )
        }

        # for GO_SS
        if(is(object,"GO_clusters")){

            # color labels of cutting tree
            colors= unique(
                get_leaves_attr(
                    slot(object,"dendrograms")$GO,
                    "edgePar"
                )
            )

            # extract the GO. cluster column
            count<-slot(slot( object,"enrich_GOs"),"data")[,"GO.cluster",with=FALSE]

            # count the number of terms by clusters
            count<-count[,list("nb"=.N),by="GO.cluster"]

            # add count to res.mds
            res.mds=merge(
                res.mds,
                count,
                by="GO.cluster",
                sort=FALSE
            )

            # for terms
            if(type=="GOterms"){

                # create trace
                p<-add_markers(
                    p,
                    data=res.mds,
                    x=~Dim.1,
                    y=~Dim.2,
                    color=~GO.cluster,
                    text = ~paste("cluster:",GO.cluster,"<br>GO.ID:",GO.ID,"<br>GO.name:",term),
                    showlegend=TRUE,
                    colors=colors,
                    marker =list(
                        size =20,
                        opacity = 0.4
                    )
                )

                # add custom layout with dropdown menu
                p<-layout(
                    p,

                    # add title
                    title=paste(measures,"distance MultiDimensional Scaling plot"),

                    # increase font size
                    font=list(size=14),

                    # add axis legends
                    xaxis=list(title="Dimension 1"),
                    yaxis=list(title="Dimension 2")
                )
            }else{

                # add count to text
                res.mds[,`:=`(text=paste(sub("<br>.+$","",text),"<br>GO.count:",res.mds$nb,sub("^.+<br>GO.ID","<br>GO.ID",text)))]

                # add trace to plot by measure
                for(x in seq_along(measures)){

                # default visualization
                if(x==1){visible=TRUE}else{visible=FALSE}

                # create trace
                p<-add_markers(
                    p,
                    data=res.mds[measures[x],on="measure"],
                    x=~Dim.1,
                    y=~Dim.2,
                    name=measures[x],
                    text=~text,
                    showlegend=FALSE,
                    sizes=c(20,50),
                    size=~nb,
                    marker =list(
                        sizemode = 'diameter',
                        opacity = 0.4,
                        line=list(color=colors),
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
          font=list(size=14),

          #################
          # add axis legends
          xaxis=list(title="Dimension 1"),
          yaxis=list(title="Dimension 2"),

          #################
          # add scrolling menu for availables measures
          updatemenus = list(

            #################
            # measure dropdown menu
            list(

              #################
              # button position
              x = 0.1,
              y = 1.1,

              buttons = lapply(seq_along(measures),function(x){

                #################
                # init visibility to all FASLE
                values=rep(FALSE,length(measures))

                #################
                # turn to true for each measures
                values[x]<-TRUE

                #################
                # output button
                list(method ="restyle",
                args = list("visible",as.list(values)),
                label=measures[x])
              })
            )
          )
        )
      }
        }

        # return or print
        if(is.null(file)){

            # return the plot
            p

        }else{

            # print heatmap
            export(p,file=file)
        }
    }
)
