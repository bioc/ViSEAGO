###################
# package startup check
.onAttach<-function(libname, pkgname) {

  ###################
  # package startup check
  toinstall=c(
    "data.table",
    "AnnotationDbi",
    "AnnotationForge",
    "biomaRt",
    "dendextend",
    "DiagrammeR",
    "GOSemSim",
    "DT",
    "dynamicTreeCut",
    "ggplot2",
    "GO.db",
    "heatmaply",
    "htmltools",
    "igraph",
    "methods",
    "org.Mm.eg.db",
    "plotly",
    "topGO",
    "Rcpp",
    "RColorBrewer",
    "Rgraphviz",
    "R.utils",
    "scales",
    "UpSetR",
    "webshot",
    "BiocStyle",
    "knitr",
    "rmarkdown",
    "corrplot"
  )

  ###################
  # previously installed
  check=toinstall%in%utils::installed.packages()[,1]

  ###################
  # install if needed
  if(!base::all(check)){
    if (!requireNamespace("BiocManager")) install.packages("BiocManager")
    BiocManager::install(toinstall[!check],dependencies=T)
    webshot::install_phantomjs()
  }
}
