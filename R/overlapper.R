#' @title build all sets combinations intersections
#' @description  build all sets combinations intersections
#' @family enrich_GO_terms
#' @family visualization
#' @details This internal function use build all sets combinations intersections needed for \code{\link{upset}}.
#' @return a \code{\link[base]{list}}.
#' @examples
#' \dontrun{
#' ###################
#' # build all intersections combinations
#' ViSEAGO::overLapper(setlist)
#' }
#' @keywords internal
#' @export
overLapper<-function(setlist){

    # Create intersect matrix (removes duplicates!)
    complexity=seq_along(setlist)
    setunion <- sort(unique(unlist(setlist)))
    setmatrix <- vapply(names(setlist), function(x){
        setunion %in% unique(setlist[[x]])},
        rep(TRUE,length(setunion))
    )
    rownames(setmatrix) <- setunion
    storage.mode(setmatrix) <- "numeric"

    # Create all possible sample combinations within requested complexity levels
    labels <- names(setlist)
    allcombl <- lapply(complexity, function(x){
        combn(labels, m=x, simplify=FALSE)
    })
    allcombl <- unlist(allcombl, recursive=FALSE)
    complevels <- vapply(allcombl, length,0)

    # vennSets function
    vennSets <- function(setmatrix=setmatrix, allcombl=allcombl, index=1) {
        mycol1 <-which(colnames(setmatrix) %in% allcombl[[index]])
        mycol2 <-which(!colnames(setmatrix) %in% allcombl[[index]])
        cond1 <-rowSums(setmatrix[, rep(mycol1, 2)]) == 2 * length(mycol1)
        cond2 <-rowSums(setmatrix[, rep(mycol2, 2)]) == 0
        return(setunion[cond1 & cond2])
    }

    # apply vennSets to all conbinations
    vennOLlist <- lapply(seq(along=allcombl), function(x){
        vennSets(setmatrix=setmatrix, allcombl=allcombl, index=x)
    })
    names(vennOLlist) <- vapply(allcombl, paste, collapse="-","")

    # returnvennOLlist
    return(vennOLlist)
}
