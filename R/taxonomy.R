#' @title Display Organism Scientific and common name from taxid.
#' @description  Display Organism Scientific and common name from taxid.
#' @importFrom data.table data.table
#' @family genomic_ressource
#' @details This internal function use \href{https://www.ncbi.nlm.nih.gov/books/NBK25500/}{E-utilities} summary to display in
#'  \code{\link[data.table]{data.table}} Organism Scientific name and common name, form a \code{\link[base]{vector}} of taxid.
#' @return a \code{\link[data.table]{data.table}}.
#' @examples
#' ###################
#' # Organism Scientific and common name from taxid
#' ViSEAGO::taxonomy(taxid="9031")
#' @keywords internal
#' @export
taxonomy=function(...){

  ###################
  # taxonomy ids
  taxid=base::unique(...)

  ###################
  # internal function for pattern  extraction
  pattern.extract=function(query,m){

    ###################
    # for each query line
    Data=base::lapply(seq_along(query),function(i){

      ###################
      # if not empty m match
      if(base::length(stats::na.omit(m[[i]][1]))>0){

        ###################
        # extract capture.start argument
        capture=base::attr(m[[i]],"capture.start")

        ###################
        # extract corresponding values in query
        capture=base::substring(query[i], capture,capture+base::attr(m[[i]],"capture.length")-1)

        ###################
        # convert in data.table
        capture=data.table::data.table(base::t(capture))

      }else{

        ###################
        # else NA
        NA
      }
    })

    ###################
    # rbind.data.table
    Data<-data.table::rbindlist(Data)

    ###################
    # add header
    base::colnames(Data)<-base::attr(m[[1]],"capture.name")

    ###################
    # replace "" or \t values by NA
    Data[CommonName%in%c("","\t"),CommonName:=NA]

    ###################
    # return query
    Data
  }

  ###################
  # core address
  core="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?version=2.0&db=taxonomy"

  ###################
  # create submission query
  query <-base::paste(core,"&id=",base::paste(taxid,collapse=","),sep = "")

  ###################
  # submit and retrieve
  query=base::paste(base::scan(query, what ="",sep="\n",quiet = T),collapse="")

  ###################
  # parse results
  query<-base::substring(query,base::unlist(gregexpr("<DocumentSummary ",query)),
  base::unlist(base::gregexpr("</DocumentSummary>",query)))

  ###################
  # extraction pattern
  pattern=c("<DocumentSummary uid=\"(?<taxid>[[:digit:]]*)\"",
  ".*<ScientificName>(?<ScientificName>.*)</ScientificName>",
  "\t<CommonName>(?<CommonName>.*)</CommonName>.*")

  ###################
  # find pattern
  m=base::gregexpr(base::paste(pattern,collapse=""),query,perl=TRUE)

  ###################
  # extract  results in data.frame and return
  pattern.extract(query,m)
}
