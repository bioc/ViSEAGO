#' @title Display Organism Scientific and common name from taxid.
#' @description  Display Organism Scientific and common name from taxid.
#' @importFrom data.table data.table rbindlist
#' @family genomic_ressource
#' @details This internal function use \href{https://www.ncbi.nlm.nih.gov/books/NBK25500/}{E-utilities} summary to display in
#'  \code{\link[data.table]{data.table}} Organism Scientific name and common name, form a \code{\link[base]{vector}} of taxid.
#' @return a \code{\link[data.table]{data.table}}.
#' @examples
#' ###################
#' # Organism Scientific and common name from taxid
#' Data<-ViSEAGO::taxonomy("9031")
#' @keywords internal
#' @export
taxonomy=function(...){

    # taxonomy ids
    taxid=unique(...)

    # internal function for pattern  extraction
    pattern.extract<-function(query,m){

        # for each query line
        Data=lapply(seq_along(query),function(i){

            # if not empty m match
            if(length(stats::na.omit(m[[i]][1]))>0){

                # extract capture.start argument
                capture=attr(m[[i]],"capture.start")

                # extract corresponding values in query
                capture=substring(
                    query[i],
                    capture,
                    capture+attr(m[[i]],"capture.length")-1
                )

                # convert in data.table
                data.table(t(capture))

            }else{

                # else NA
                NA
            }
        })

        # rbind.data.table
        Data<-rbindlist(Data)

        # add header
        colnames(Data)<-attr(m[[1]],"capture.name")

        # replace "" or \t values by NA
         Data[c("","\t"),"CommonName":="NA",on="CommonName"]

        # return query
        return(Data)
    }

    # core address
    core="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?version=2.0&db=taxonomy"

    # create submission query
    query <-paste(
        core,
        "&id=",
        paste(taxid,collapse=","),
        sep = ""
    )

  # submit and retrieve
    query=paste(
        scan(
            query,
            what ="",
            sep="\n",
            quiet = TRUE
        ),
        collapse=""
    )

    # parse results
    query<-substring(
        query,
        unlist(gregexpr("<DocumentSummary ",query)),
        unlist(gregexpr("</DocumentSummary>",query))
    )

    # extraction pattern
    pattern=c("<DocumentSummary uid=\"(?<taxid>[[:digit:]]*)\"",
    ".*<ScientificName>(?<ScientificName>.*)</ScientificName>",
    "\t<CommonName>(?<CommonName>.*)</CommonName>.*")

    # find pattern
    m=gregexpr(
        paste(pattern,collapse=""),
        query,
        perl=TRUE
    )

    # extract  results in data.frame and return
    pattern.extract(query,m)
}
