##
#' Create a peptide collection
#'
#' Constructor to create peptide collection such as the datasets available in PEP.db.
#'
#' @param ranges An \code{RangedData} object. The object should have a peptide column.
#' 
#' @seealso \code{\link{RangedData}}
#' 
#' @examples
#' 
#'
#' @importFrom IRanges RangedData DataFrame IRanges
#' @importFrom pepStat makeZpepMatrix
##

create_db <- function(rd){
  if(class(rd) != "RangedData"){
    stop("`rd' must be an object of class 'RangedData'. Given: ",class(rd))
  }
  if(rd$peptu
  zs <- makeZpepMatrix(as.character(db$peptide))
  rownames(zs) <- NULL
  values(db)[[1]] <- cbind(values(db)[[1]], DataFrame(zs))
  rownames(db) <- as.character(db$peptide)
  return(db)
:w
}
