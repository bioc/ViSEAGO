#' @title Display the ViSEAGO package diagram.
#' @description This function build the ViSEAGO package diagram available displayed in the vignette.
#' @importFrom DiagrammeR grViz
#' @param x rectangles colors for the five subgraphs with all white by defaults
#' @details This function use \pkg{DiagrammeR} package \code{\link[DiagrammeR]{grViz}} to display the ViSEAGO methods diagram with colored or not subgraphs rectangles.
#' @references
#' Knut Sveidqvist, Mike Bostock, Chris Pettitt, Mike Daines, Andrei Kashcha and Richard Iannone (2017). DiagrammeR: Create Graph Diagrams and Flowcharts Using R.
#' R package version 0.9.0. https://CRAN.R-project.org/package=DiagrammeR.
#' @return an html diagram
#' @examples
#' \dontrun{
#' ###################
#' # print diagram without focus
#' ViSEAGO::pkgdiagram(x=c("white","white","white","white","white"))
#'
#' ###################
#' # print diagram with focus one the gene list
#' ViSEAGO::pkgdiagram(x=c("black","white","white","white","white"))
#'
#' ###################
#' # print diagram with a focus on the annotation step
#' ViSEAGO::pkgdiagram(x=c("white","black","white","white","white"))
#'
#' ###################
#' # print diagram with a focus on the enrichment step
#' ViSEAGO::pkgdiagram(x=c("white","white","black","white","white"))
#'
#' ###################
#' # print diagram with a focus on the Semantic Similarity step
#' ViSEAGO::pkgdiagram(x=c("white","white","white","black","white"))
#'
#' ###################
#'  # print diagram with focus on the visualization step
#' ViSEAGO::pkgdiagram(x=c("white","white","white","white","black"))
#' }
#' @keywords internal
#' @export
pkgdiagram<-function(x=c("white","white","white","white","white")){

  ###################
  # print the diagram
  DiagrammeR::grViz(diagram=base::paste("

    ###################
    # init the graph
    digraph ViSEAGO_map {

    ###################
    # a 'graph' statement
    graph [layout = dot,overlap = false, fontsize = 10, splines= true, shape =oval,fixedsize = true,width = 2.5,fontname = Helvetica]

    subgraph cluster0 {
      label='List(s) of genes'
      fontsize = 25
      fontcolor = MistyRose2
      labelloc=t
      color=",x[1],"
      node [fillcolor=MistyRose2,style=filled]
      'Genes of interest \n and background'
    }

    node [fillcolor=orange,style=filled]
    subgraph cluster1 {
      label='Genomic ressources'
      fontsize = 25
      fontcolor = orange
      labelloc=t
      color=",x[2],"
      Bioconductor2GO
      EntrezGene2GO
      Ensembl2GO
      Uniprot2GO
      available_organisms
      annotate
    }

    node [fillcolor=lightblue1,style=filled]
    subgraph cluster2 {
      label='Enrichment tests'
      fontsize = 25
      fontcolor =lightblue1
      labelloc=t
      color=",x[3],"
      create_topGOdata
      RunTest
      merge_enrich_terms
    }

    node [fillcolor=LightCoral,style=filled]
    subgraph cluster3 {
      label='Visualization'
      fontsize = 25
      fontcolor = LightCoral
      labelloc=t
      color=",x[5],"
      show_table
      show_heatmap
      GOcount
      Upset
      MDSplot
      GOterms_heatmap
      GOclusters_heatmap
      compare_clusters
    }

    node [fillcolor=LimeGreen,style=filled]
    subgraph cluster4 {
      label='GO Semantic Similarities'
      fontsize = 25
      fontcolor = LimeGreen
      labelloc=t
      color=",x[4],"
      build_GO_SS
      compute_SS_distances
    }

    ###################
    #  'edge' statements
    'Genes of interest \n and background'-> create_topGOdata
    Bioconductor2GO->annotate
    EntrezGene2GO->annotate
    Ensembl2GO->annotate
    Uniprot2GO->annotate
    annotate->create_topGOdata
    create_topGOdata-> {RunTest merge_enrich_terms}
    RunTest->merge_enrich_terms
    merge_enrich_terms->{show_table GOcount Upset}
    merge_enrich_terms->build_GO_SS
    build_GO_SS->compute_SS_distances
    compute_SS_distances->{MDSplot GOterms_heatmap GOclusters_heatmap}
    GOterms_heatmap->{compute_SS_distances show_table show_heatmap compare_clusters}
    GOclusters_heatmap->{show_table show_heatmap MDSplot}
    Bioconductor2GO->available_organisms
    EntrezGene2GO->available_organisms
    Ensembl2GO->available_organisms
    Uniprot2GO->available_organisms

    }",sep="")
  )
}
