#' @title Build a clustering heatmap on GO terms.
#' @description This method computes a clustering heatmap based on GO terms semantic similarity.
#' @importFrom data.table data.table
#' @importFrom graphics text
#' @importFrom stats as.dist end start
#' @importFrom methods setGeneric setMethod new slot is signature
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom plotly layout
#' @importFrom heatmaply heatmaply
#' @importFrom dendextend branches_attr_by_clusters  get_nodes_attr set
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales rescale
#' @family GO_terms
#' @family semantic_similiarity
#' @family visualization
#' @param myGOs a \code{\link{GO_SS-class}} object from \code{\link{compute_SS_distances}}.
#' @param showIC \code{logical} (default to TRUE) to display the GO terms Information Content (IC) side bar.
#' @param showGOlabels \code{logical} (default to TRUE) to display the GO terms ticks on y axis.
#' @param GO.tree a named \code{list} of parameters to build and cut the GO terms \code{dendrogram}.
#'  \describe{
#'      \item{tree (a named \code{list} with:)}{
#'          \describe{
#'              \item{distance ("Wang" by default)}{distance computed from the semantic similarity which could be
#'               IC-based ("Resnik", "Rel", "Lin", or "Jiang") or graph-based ("Wang").
#'              }
#'               \item{aggreg.method ("ward.D2" by default)}{aggregation method criteria from \code{\link[stats]{hclust}} ("ward.D",
#'               "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid") to build a \code{dendrogram}.}
#'              \item{rotate}{sort the branches of the tree based on a vector - eithor of labels order or the labels in their new order}
#'          }
#'      }
#'      \item{cut (a named \code{list} with:)}{
#'          \describe{
#'              \item{static (default to NULL)}{a \code{numeric} value that is the height (between 0 and 1),
#'              or the number of clusters (value > 1) to cut the \code{dendrogram}.
#'              }
#'               \item{dynamic (a named \code{list} which only contains \code{\link[dynamicTreeCut]{cutreeDynamic}}
#'               options values below)
#'               }{
#'                 \describe{
#'                     \item{pamStage (default to TRUE)}{second (PAM-like) stage will be performed.}
#'                     \item{pamRespectsDendro (default to TRUE)}{PAM stage will respect the dendrogram in the sense that objects
#'                      and small clusters will only be assigned to clusters that belong to the same branch that the objects or small
#'                      clusters being assigned belong to.}
#'                     \item{deepSplit (default to 2)}{provides a rough control over sensitivity for cluster splitting (range 0 to 4).
#'                       The higher the value (or if TRUE), the more and smaller clusters will be produced.}
#'                     \item{minClusterSize (default to 2)}{minimum cluster size.}
#'                 }
#'              }
#'          }
#'      }
#' }
#' @param samples.tree a named \code{list} of parameters to build and cut the samples \code{dendrogram} (default to NULL).
#'  \describe{
#'      \item{tree (a named \code{list} with:)}{
#'          \describe{
#'              \item{distance ("pearson" by default)}{distance computed that could be correlation ("abs.pearson","pearson", "kendall", or "spearman"),
#'              or \code{dist} method (euclidean", "maximum", "manhattan", "canberra", "binary", or "minkowski).}
#'               \item{aggreg.method ("average" by default)}{same options than for \code{GO.tree} argument}
#'          }
#'      }
#'      \item{cut}{same options than for \code{GO.tree} argument.}
#' }
#' @details This method computes a clustering heatmap based on GO terms semantic similarity (computed with \code{\link{compute_SS_distances}}).\cr
#' The dendrogram produced could be cutted in static or dynamic mode.\cr
#'  \enumerate{
#'   \item build dendrograms on  GO terms and optionally on samples.
#'   \item cut in static or dynamic mode and color the dendrogram branchs.
#'   \item build an interactive clustering heatmap based on \code{\link[heatmaply]{heatmaply}}.
#'  }
#' @return a \code{\link{GO_clusters-class}} object.
#' @references
#' Matt Dowle and Arun Srinivasan (2017). data.table: Extension of `data.frame`. R package version 1.10.4. https://CRAN.R-project.org/package=data.table.
#'
#' Tal Galili (2015). dendextend: an R package for visualizing, adjusting, and comparing trees of hierarchical clustering.
#' Bioinformatics. DOI:10.1093/bioinformatics/btv428.
#'
#' Tal Galili (2017). heatmaply: Interactive Cluster Heat Maps Using 'plotly'.
#' R package version 0.9.1. https://CRAN.R-project.org/package=heatmaply.
#'
#' Peter Langfelder, Bin Zhang and with contributions from Steve Horvath (2016). dynamicTreeCut: Methods for Detection of Clusters
#' in Hierarchical Clustering Dendrograms. R package version 1.63-1. https://CRAN.R-project.org/package=dynamicTreeCut.
#'
#' Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1-2. https://CRAN.R-project.org/package=RColorBrewer.
#'
#' Carson Sievert, Chris Parmer, Toby Hocking, Scott Chamberlain, Karthik Ram, Marianne Corvellec and Pedro Despouy (2017).
#' plotly: Create Interactive Web Graphics via 'plotly.js'. R package version 4.6.0. https://CRAN.R-project.org/package=plotly.
#'
#' Hadley Wickham (2016). scales: Scale Functions for Visualization. R package version 0.4.1. https://CRAN.R-project.org/package=scales.
#'
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
#' @include GO_clusters.R
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
#'  distance="Wang"
#' )
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
#' }
#' @name GOterms_heatmap
#' @rdname GOterms_heatmap-methods
#' @exportMethod GOterms_heatmap
setGeneric(name="GOterms_heatmap",def=function(
  myGOs,
  showIC=TRUE,
  showGOlabels=TRUE,
  GO.tree=base::list(
    tree=base::list(
      distance="Wang",
      aggreg.method="ward.D2",
      rotate=NULL),
    cut=base::list(
      dynamic=base::list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE
        ,deepSplit=2,
        minClusterSize =2
      )
    )
  ),
  samples.tree=NULL
){
  base::standardGeneric("GOterms_heatmap")
})

#' @rdname GOterms_heatmap-methods
#' @aliases GOterms_heatmap
setMethod("GOterms_heatmap",
  methods::signature(
    myGOs="GO_SS"
  ),
  definition=function(myGOs,showIC,showGOlabels,GO.tree,samples.tree){

  ###################
  # check entry
  ###################
  if(base::length(methods::slot(myGOs,"terms_dist"))==0){
    base::stop("Please compute Semantic Similarity distance with ViSEAGO::compute_SS_distances()")
  }

  ###################
  # creates trees
  ###################

    ###################
    # add default params if not available in the tree list
    Tree.params<-function(...){

      ###################
      # keep the name
      name<-base::as.character(...)

      ###################
      # Tree
      Tree<-list(base::get(name))[[1]]

      if(!base::is.null(Tree)){

        ###################
        # extract list names
        tree<-base::names(Tree$tree)

        ###################
        # default distance if not available
        if("distance"%in%tree){
          distance<-Tree$tree$distance
        }else{
          if(name=="GO.tree"){
            distance="Wang"
          }else{
            distance="pearson"
          }
        }

        ###################
        # default distance if not available
        if("aggreg.method"%in%tree){
          aggreg.method<-Tree$tree$aggreg.method
        }else{
          if(name=="GO.tree"){
            aggreg.method="ward.D2"
          }else{
            aggreg.method="average"
          }
        }

        ###################
        # default order
        if("rotate"%in%tree){
          rotate<-Tree$tree$rotate
        }else{
          rotate=NULL
        }

        ###################
        # extract list names
        if(!base::is.null(Tree$cut)){

          ###################
          # default static cut tree method if not available
          tree=base::names(Tree$cut)

          ###################
          # default static cut tree method if not available
          if("static"%in%tree){
            static<-Tree$cut$static
          }else{
            static=NULL
          }

          if(base::is.null(static)){

            ###################
            # extract list names
            tree<-base::names(Tree$cut$dynamic)

            ###################
            # default aggreg.method if not available
            if("pamStage"%in%tree){
              pamStage<-Tree$cut$dynamic$pamStage
            }else{
              pamStage=TRUE
            }

            ###################
            # default aggreg.method if not available
            if("pamRespectsDendro"%in%tree){
              pamRespectsDendro<-Tree$cut$dynamic$pamRespectsDendro
            }else{
              pamRespectsDendro=TRUE
            }

            ###################
            # default aggreg.method if not available
            if("deepSplit"%in%tree){
              deepSplit<-Tree$cut$dynamic$deepSplit
            }else{
              deepSplit=2
            }

            ###################
            # default aggreg.method if not available
            if("minClusterSize"%in%tree){
              minClusterSize<-Tree$cut$dynamic$minClusterSize
            }else{
                minClusterSize =2
            }

            ###################
            # build dynamic
            dynamic=base::list(
              pamStage=pamStage,
              pamRespectsDendro=pamRespectsDendro,
              deepSplit=deepSplit,
              minClusterSize =minClusterSize
            )
          }else{

            ###################
            # build dynamic
            dynamic=NULL
          }


          ###################
          # build cut list
          cut=base::list(
            static=static,
            dynamic=dynamic
          )
        }else{

          ###################
          # build empty cut list
          cut=NULL
        }

        ###################
        # Tree with default params if nessecary
        base::list(
          tree=base::list(
            distance=distance,
            aggreg.method=aggreg.method,
            rotate=rotate
          ),
          cut=cut
        )
      }else{
        NULL
      }
    }

    ###################
    # row.tree with default if NA
    row.tree<-Tree.params("GO.tree")

    ###################
    # col.tree with default if NA
    col.tree<-Tree.params("samples.tree")

    ###################
    # check dendrogram entry
    if(base::is.null(row.tree$tree$distance)){
      base::stop(paste("please enter a myGOs object computed SS distance name (",
      paste(
        base::names(
          methods::slot(
            myGOs,
            "terms_dist"
            )
          ),
        collapse=", "),
      ") in the GO.tree distance argument",
      sep="")
      )
    }

  ###################
  # load Data
  ###################

    ###################
    # data input
    sResults<-methods::slot(
      methods::slot(
        myGOs,
        "enrich_GOs"),
      "data"
    )

    ###################
    # extract pvalues columns
    mat<-base::as.matrix(sResults[,grep("-log10.pvalue",names(sResults)),with=FALSE])

    ###################
    # remove -log10.pvalue in colnames
    base::colnames(mat)<-base::gsub("\\.-log10.pvalue","",base::colnames(mat))

    ###################
    # duplicate sResults in Data
    Data<-sResults

    ###################
    #  add modfied terms definition
    base::row.names(mat)<-base::paste("<br>GO.ID:",Data$GO.ID,"<br>GO.name:",Data$term)

    ###################
    # if Inf in -log10pvalue, set 9
    for(i in  base::seq_len(base::ncol(mat))){mat[mat[,i]=="Inf",i]<-9}

    ###################
    #  remove Data
    base::rm(Data)

    ###################
    # transpose for GO in rows and samples in columns  in the levelplot
    x<-base::t(mat)

    ###################
    # Add IC in panes
    if(showIC==TRUE){

      ###################
      # extract values
      IC<-base::round(methods::slot(myGOs,"IC")[sResults$GO.ID],digits=2)
    }

  ###################
  # distance correlation
  ###################

    ###################
    # distance function
    dist.fun<-function(...){

      ###################
      # text tree
      Tree<-base::as.character(...)

      ###################
      # row or col tree
      if(base::length(base::grep("col",Tree))==1){col=T}else{col=F}

      ###################
      # get tree
      Tree<-base::get(Tree)

      ###################
      # if tree on column
      if(!base::is.null(Tree)){

        ###################
        # column tree with correlation distance
        if(base::is.character(Tree$tree$distance) || Tree$tree$distance%in%c("abs.pearson","pearson", "kendall","spearman")){

          ###################
          # if tree on column
          if(col==T){x<-base::t(x)}

          ###################
          # for absolute pearson correlation
          if(Tree$tree$distance=="abs.pearson"){

            ###################
            # col.tree with absolute pearson correlation
            stats::as.dist(1-base::abs(stats::cor(x,method="pearson")))

          }else{

            ###################
            # other correlation distances
            stats::as.dist(1-stats::cor(x,method=Tree$tree$distance))
          }

        ###################
        # column tree without correlation distance
        }else{

          ###################
          # if row on column
          if(col==F){x<-t(x)}

          ###################
          # calculate distance
          stats::dist(x,method=Tree$tree$distance)
        }
      }
    }

    ###################
    # columns distance
    col.dist<-dist.fun("col.tree")

    ###################
    # row distance (precalculated distance based on similarity of GO terms)
    row.dist<-methods::slot(myGOs,"terms_dist")[[row.tree$tree$distance]]

  ###################
  # aggregation
  ###################

    ###################
    # aggregation function
    aggreg.fun<-function(...){

      ###################
      # text tree
      Tree<-base::as.character(...)

      ###################
      # add calculated distance
      dist<-base::get(base::paste(base::substring(Tree,1,3),"dist",sep="."))

      ###################
      # if distance matrix not empty
      if(!base::is.null(dist)){

        ###################
        # get tree
        Tree<-base::get(Tree)

        ###################
        # perform hclust
        hc<-stats::hclust(dist, method =Tree$tree$aggreg.method)

        ###################
        # rotate
        if(!base::is.null(Tree$tree$rotate)){
          hc<-dendextend::rotate(hc,Tree$tree$rotate)
        }

        ###################
        # return roated hc
        base::return(hc)
      }
    }

    ###################
    # columns distance
    col_tree<-aggreg.fun("col.tree")

    ###################
    # row distance
    row_tree<-aggreg.fun("row.tree")

  ###################
  # create dendrogram(s)
  ###################

    ###################
    # visible binding for global variable
    # Bioconductor submission "checking R code possible problems" step
    col.ord<-NULL
    row.ord<-NULL

    ###################
    # create row and column dendrograms
    for(i in c("row_tree","col_tree")){

      ###################
      # create dendrogram
      if(!base::is.null(base::get(i))){

        ###################
        # dendrogram name
        dd<-base::paste("dd",base::substring(i,1,3),sep=".")

        ###################
        # dendrogram convert
        dend<-stats::as.dendrogram(base::get(i))

        ###################
        # create dendrogram
        base::assign(dd,dend)

        ###################
        # create ordering vector
        base::assign(
          base::paste(
            base::substring(i,1,3),
            "ord",
            sep="."
          ),
          stats::order.dendrogram(
            base::get(dd)
          )
        )

     }else{

        ###################
        # ordering vector ni null according column and row
        if(i=="col_tree"){Dat<-base::seq_len(base::nrow(x))}else{Dat<-base::seq_len(base::ncol(x))}

        ###################
        # create ordering vector
        base::assign(base::paste(base::substring(i,1,3),"ord",sep="."),Dat)
      }
    }

  ###################
  # cut trees
  ###################

    ###################
    # cut tree
    cut.tree<-function(...){

      ###################
      # text tree
      Tree<-base::as.character(...)

      ###################
      # text tree
      if(!base::is.null(base::get(Tree))){

        ###################
        # add calculated distance
        dendro=base::get(sub("\\.","_",Tree))

        ###################
        # add calculated distance
        dist<-base::get(base::paste(base::substring(Tree,1,3),"dist",sep="."))

        ###################
        # add order
        ord<-base::get(base::paste(base::substring(Tree,1,3),"ord",sep="."))

        ###################
        # get tree
        Tree<-base::get(Tree)

          ###################
          # cut or not
          if(!base::is.null(Tree$cut)){

            ###################
            # cut dynamic
            if(base::is.null(Tree$cut$static)){

              ###################
              # cut dynamic
              gp<-dynamicTreeCut::cutreeDynamic(
                dendro =dendro,
                method="hybrid",
                verbose=F,
                distM =as.matrix(dist),
                pamStage=Tree$cut$dynamic$pamStage,
                pamRespectsDendro=Tree$cut$dynamic$pamRespectsDendro,
                deepSplit=Tree$cut$dynamic$deepSplit,
                minClusterSize =Tree$cut$dynamic$minClusterSize
              )

              ###################
              # groups
              clust=unique(gp[ord])

              ###################
              # new clusters names table
              clust=data.table::data.table(ini=clust,new=base::seq_len(length(clust)))

              ###################
              # convert gp to data.table
              gp<-data.table::data.table(ini=gp)

              ###################
              # merge gp with new clusters names
              gp<-merge(gp,clust,by="ini",all.x=T,sort=F)

              ###################
              # extract only new clusters names in the initial order
              gp<-gp$new

              ###################
              # add go names
              names(gp)<-base::attr(dist,"Labels")

              ###################
              # return the object
              gp

          }else{

            ###################
            # static cut tree by group numbers
            if(Tree$cut$static>1){

              ###################
              # cut row tree
              stats::cutree(dendro, k =Tree$cut$static,h = NULL)[ord]

            }else{

              ###################
              # cut row tree
              stats::cutree(dendro,k=NULL,h=Tree$cut$static)[ord]
            }
          }
        }else{

          ###################
          # text tree
          Tree<-base::as.character(...)

          ###################
          # unique group
          gp<-base::rep(1,base::length(base::get(base::sub("tree","ord",Tree))))

          ###################
          # give col names attributes
          if(Tree=="col.tree"){base::attr(gp,"names")<-base::row.names(x)}

          ###################
          # give row names attributes
          if(Tree=="row.tree"){base::attr(gp,"names")<-base::colnames(x)}

          ###################
          # return gp
          gp
        }

      }else{

        ###################
        # unique group
        gp<-base::rep(1,length(base::get(base::sub("tree","ord",Tree))))

        ###################
        # give col names attributes
        if(Tree=="col.tree"){base::attr(gp,"names")<-base::row.names(x)}

        ###################
        # give row names attributes
        if(Tree=="row.tree"){base::attr(gp,"names")<-base::colnames(x)}

        ###################
        # return gp
        gp
      }
    }

    ###################
    # cut row tree
    row.gp<-cut.tree("row.tree")

    ###################
    # cut column tree
    col.gp<-cut.tree("col.tree")

    ###################
    # add cluster to row.names
    base::row.names(mat)<-base::paste("<br>cluster:",row.gp,base::row.names(mat))

  ###################
  # color row dendrogram according clusters
  ###################

    ###################
    # color labels of cutting tree
    palette=base::c("darksalmon","cyan","orchid","midnightblue","hotpink",
      "sienna","aquamarine","royalblue","forestgreen","salmon","maroon",
      "green","deepskyblue","olivedrab","springgreen","limegreen","magenta",
      "peru","tan","steelblue","orange","burlywood","skyblue","red",
      "darkslateblue","plum","goldenrod","cadetblue","indianred","pink",
      "slateblue","turquoise","darkorchid","cornflowerblue","darkolivegreen",
      "orangered","powderblue","rosybrown","darkgoldenrod","seagreen",
      "darkturquoise","thistle","dodgerblue","sandybrown","purple","tomato",
      "deeppink", "darkmagenta","darkred","darkkhaki","coral","lawngreen",
      "darkorange","black","chocolate","firebrick","darkseagreen","blue",
      "navy","darkgreen"
    )

    ###################
    # color branches according clusters
    dd.row<-dendextend::branches_attr_by_clusters(dd.row,row.gp[row.ord], values =  palette)

    ###################
    # color labels of cutting tree
    colors=stats::na.omit(base::unique(base::unlist(dendextend:: get_nodes_attr(dd.row,"edgePar"))))

    ###################
    # colors table
    colors=data.table::data.table(gp=base::unique(row.gp[row.ord]),color=colors)

    ###################
    # merge cluster term assignation and corresponding color
    colors<-merge(data.table::data.table(gp=row.gp[row.ord]),colors,by="gp",all.x=T,sort=F)

    ###################
    # assign text color
    dd.row<-dendextend::set(dd.row,"labels_col",colors$color)

    ###################
    # create dendrogram
    dd.row<-dendextend::set(dd.row,"labels_cex",1)

  ###################
  # column dendrogram according clusters
  ###################

    if(!is.null(col.dist)){

      ###################
      # color branches according clusters
      dd.col<-dendextend::branches_attr_by_clusters(dd.col,col.gp[col.ord], values =palette)

      ###################
      # color labels of cutting tree
      colors=stats::na.omit(base::unique(base::unlist(dendextend:: get_nodes_attr(dd.col,"edgePar"))))

      ###################
      # colors table
      colors=data.table::data.table(gp=base::unique(col.gp[col.ord]),color=colors)

      ###################
      # merge cluster term assignation and corresponding color
      colors<-merge(data.table::data.table(gp=col.gp[col.ord]),colors,by="gp",all.x=T,sort=F)

      ###################
      # assign text color
      dd.col<-dendextend::set(dd.col,"labels_col",colors$color)

      ###################
      # create dendrogram
      dd.col<-dendextend::set(dd.col,"labels_cex",1)
    }

  ###################
  # draw heatmap
  ###################

    ###################
    # If gene background not the same
    if(methods::slot(methods::slot(myGOs,"enrich_GOs"),"same_genes_background")==FALSE){

      ###################
      # show warning
      base::warning(
        call. =FALSE,
        "Not equal genes background in all conditions:\n--> pvalues converted to significant (1) or not (0) by condition in the cluster-heatmap"
      )

      ###################
      # reduce to significant or not (p<0.01, -log10(p)>2)
      mat<-base::apply(mat,base::c(1,2),function(x){
        if(base::is.na(x)){
          NA
        }else{
          if(x<2){0}else{1}
        }
      })
    }

    ###################
    # draw heatmapply
    hm<-heatmaply::heatmaply(

      ###################
      # the initial matrix
      x=mat,

      ###################
      # row labels
      labRow=base::row.names(mat),

      ###################
      # columns labels
      labCol=base::colnames(mat),

      ###################
      # the row dendrogram
      Rowv=dd.row,

      ###################
      # the ordered matrix according dendrograms for columns
      Colv=if(!is.null(col.dist)){dd.col}else{FALSE},

      ###################
      # the IC information
      RowSideColors=if(showIC==TRUE){data.table::data.table(IC=IC)}else{NULL},

      ###################
      # the IC information
      row_side_palette =if(showIC==TRUE){grDevices::colorRampPalette(c("#FFFFFF","#49006A"))}else{NULL},

      ###################
      # the color palette
      scale_fill_gradient_fun=if(methods::slot(methods::slot(myGOs,"enrich_GOs"),"same_genes_background")==T){

        ###################
        # if same gene background
        ggplot2::scale_fill_gradient2(
          name="-log10pvalue",
          low="white",
          mid="white",
          high ="#99000D",
          midpoint = 1.3
        )
      }else{

        ###################
        # if not same gene background
        ggplot2::scale_fill_gradient2(
          name="-log10pvalue",
          low="white",
          mid="white",
          high ="#99000D",
          midpoint =0
        )
      },

      ###################
      # the width of dendrogramm
      branches_lwd = 0.4,

      ###################
      # color bar length
      colorbar_len=0.05
    )

    ###################
    # with column dendrogram and showIC
    col=base::vapply(hm$x$data,function(x){base::length(x$text)},0)
    col<-base::which(col==(base::nrow(mat)*base::ncol(mat)))

    ###################
    # for each column
    text<-hm$x$data[[col]]$text

    ###################
    # for each column
    for (i in base::seq_len(base::ncol(text))){
      text[,i]<-base::gsub("(<br>|^)row: ","",text[,i])
    }

    ###################
    # return
    hm$x$data[[col]]$text<-text

    ###################
    # modify row_side_color hover text
    if(showIC==TRUE){

      ###################
      # if column dendrogram
      start=col+2

      ###################
      # modify row_side_color hover text
      for(i in start:(start+base::nrow(mat)-1)){

        ###################
        # IC color font correction for maxima value rounded to 0 digits
        if(base::round(base::as.numeric(hm$x$data[[i]]$name),digits=0)==base::round(base::max(IC),digits=0)){

          ###################
          # correct maxima value color purple "#49006A" in rgba
          hm$x$data[[i]]$fillcolor<-"rgba(73,0,106,1)"
        }

        ###################
        # store the text
        text<-hm$x$data[[i]]$text

        ###################
        # extract the corresponding row names value
        text<-base::row.names(mat)[
        base::as.numeric(base::gsub("^.+ ","",text))]

        ###################
        # return values
        hm$x$data[[i]]$text<-paste(
        "value:",base::round(IC[base::gsub("^.+GO.ID: | <br>GO.name.+$","",text)],digits=2),
        "<br>column: IC<br>row:",text)
      }
    }

    ###################
    # custom row text
    row.text=base::gsub("^.+GO.name: ","",
    base::rev(base::row.names(mat)[row.ord]))

    ###################
    # cut very long definition
    row.text[nchar(row.text)>50]<-base::paste(base::substring(
    row.text[nchar(row.text)>50],1,50),"...",sep="")

    ###################
    # modify layout
    hm<-plotly::layout(
      hm,

      #################
      # add title
      title=paste(row.tree$tree$distance,"GOterms distance clustering heatmap plot"),

      #################
      # title size
      font=base::list(size=14),

      #################
      # set margin
      margin =base::list(l =300,r=0, b =150,t=50)
    )

    ###################
    # modify layout
    if(!is.null(col.dist)){

      ###################
      # modify layout
      hm<-plotly::layout(
        hm,

        #################
        # domain
        yaxis=base::list(domain=base::c(0.95,1)),

        #################
        # y axis
        yaxis2=base::list(
          domain=base::c(0,0.95),
          family="Times New Roman",
          tickmode="array",
          tickvals=base::seq_len(base::nrow(mat)),
          ticktext=row.text,
          tickfont=base::list(size=10),
          showticklabels=showGOlabels
        )
      )
    }else{

      ###################
      # modify layout
      hm<-plotly::layout(
        hm,

        #################
        # y axis
        yaxis=base::list(
          family="Times New Roman",
          tickmode="array",
          tickvals=if(showGOlabels==TRUE){base::seq_len(base::nrow(mat))}else{NULL},
          ticktext=row.text,
          tickfont=base::list(size=10),
          showticklabels=showGOlabels
        )
      )
    }

    ###################
    # modify layout
    if(showIC==TRUE){

      hm<-plotly::layout(hm,

        #################
        # x axis
        xaxis=base::list(
          domain=c(0.025,0.50),
          family="Times New Roman",
          tickfont=base::list(size=10)
        ),
        xaxis2=base::list(
          domain=c(0.50,0.55)
        ),
        xaxis3=base::list(
          domain=c(0.55,1)
        )
      )
    }else{

      hm<-plotly::layout(hm,

        #################
        # x axis
        xaxis=base::list(
          domain=c(0.025,0.55),
          family="Times New Roman",
          tickfont=base::list(size=10)
        ),
        xaxis2=base::list(
          domain=c(0.55,1)
        )
      )
    }

    ###################
    #  hm to list
    hm<-base::list(hm)

    ###################
    # give names to hm list
    base::names(hm)<-"GOterms"

    ###################
    # create an empty column dendrogram if needed
    if(!"dd.col"%in%base::ls()) dd.col=NULL

  ###################
  # export results
  ###################

    ###################
    # bind data
    sResults<-data.table::data.table(GO.cluster=row.gp,IC=IC,sResults)

    ###################
    # ordering results
    sResults<- sResults[row.ord]

    ###################
    # replace data by sResults in myGOs
    methods::slot(methods::slot(myGOs,"enrich_GOs"),"data")<-sResults

    ###################
    # remove not used elements in row.tree
    if(!is.null(row.tree$cut$static)){

      ###################
      # keep only static parameters
      row.tree$cut$dynamic<-NULL

      ###################
      # remove not used elements in row.tree
      if(row.tree$cut$static<=1){

        ###################
        # add h to static
        base::names(row.tree$cut)<-base::paste(names(row.tree$cut),"(h)")

      }else{

        ###################
        # add k to static
        base::names(row.tree$cut)<-base::paste(names(row.tree$cut),"(k)")
      }
    }

    ###################
    # return heatmap
    methods::new("GO_clusters",
      db=methods::slot(myGOs,"db"),
      stamp =methods::slot(myGOs,"stamp"),
      organism=methods::slot(myGOs,"organism"),
      ont=methods::slot(myGOs,"ont"),
      topGO=methods::slot(myGOs,"topGO"),
      IC=methods::slot(myGOs,"IC"),
      enrich_GOs=methods::slot(myGOs,"enrich_GOs"),
      terms_dist=methods::slot(myGOs,"terms_dist")[row.tree$tree$distance],
      hcl_params=base::list(
        GO.tree=row.tree,
        sample.tree=col.tree
      ),
      dendrograms=base::list(samples=dd.col,GO=dd.row),
      samples.gp=col.gp,
      heatmap=hm
    )
})
