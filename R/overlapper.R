#' @title  build all sets combinations intersections
#' @description   build all sets combinations intersections
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
overLapper <- function(setlist) {

  ###################
  # Create intersect matrix (removes duplicates!)
  complexity=base::seq_len(base::length(setlist))
  setunion <- base::sort(base::unique(base::unlist(setlist)))
  setmatrix <- base::vapply(base::names(setlist), function(x){
    setunion %in% base::unique(setlist[[x]])},
    base::rep(TRUE,base::length(setunion))
  )
  base::rownames(setmatrix) <- setunion
  base::storage.mode(setmatrix) <- "numeric"

  ###################
  # Create all possible sample combinations within requested complexity levels
  labels <- base::names(setlist)
  allcombl <- base::lapply(complexity, function(x) utils::combn(labels, m=x, simplify=FALSE))
  allcombl <- base::unlist(allcombl, recursive=FALSE)
  complevels <- base::vapply(allcombl, length,0)

  ###################
  # vennSets function
  vennSets <- function(setmatrix=setmatrix, allcombl=allcombl, index=1) {
    mycol1 <-base::which(base::colnames(setmatrix) %in% allcombl[[index]])
    mycol2 <-base::which(!base::colnames(setmatrix) %in% allcombl[[index]])
    cond1 <-base::rowSums(setmatrix[, base::rep(mycol1, 2)]) == 2 * base::length(mycol1)
    cond2 <-base::rowSums(setmatrix[, base::rep(mycol2, 2)]) == 0
    return(setunion[cond1 & cond2])
  }

  ###################
  # apply vennSets to all conbinations
  vennOLlist <- base::lapply(seq(along=allcombl), function(x){
    vennSets(setmatrix=setmatrix, allcombl=allcombl, index=x)
  })
  base::names(vennOLlist) <- vapply(allcombl, paste, collapse="-","")

  ###################
  # returnvennOLlist
  base::return(vennOLlist)
}
