#' @name pep_hxb2
#' @title Reference peptide dataset for HIV enveloppe
#'
#' @description A peptide dataset created from the multiple alignment of the 
#' reference sequence hxb2 and the seven subtypes A, B, C, D, M, CRF01, CRF02
#' 
#' @usage data(pep_hxb2)
#' @format A \code{GRanges} object. One row per peptide.
#'
#' @details
#' \itemize{
#'   \item{peptides:}{ The rownames are 15mers from the envelope of 7 subtypes of HIV-1.}
#'   \item{seqnames:}{ An information regarding the genomic location of the peptide.}
#'   \item{ranges:}{ The start and end coordinate of the peptide based on the coordinates of the reference hxb2.}
#'   \item{aligned:}{ The aligned column represents the peptide as it is aligned in the multiple alignment.}
#'   \item{trimmed:}{ The trimmed column represents the part of hxb2's sequence aligned with the peptide.}
#'   \item{zSums:}{ The sum of the zScale of each amino-acid of the 15mer.}
#'   \item{clade:}{ A list of the subtypes which contain this peptide.}
#' }
#'
#' @docType data
#' @importClassesFrom GenomicRanges GRanges
#' @note The alignments used to create this collections are available in the 
#' extdata/alignments folder, as well as the raw sequences on fasta format. 
#' @references http://www.hiv.lanl.gov/content/sequence/HIV/REVIEWS/HXB2.html
##
NULL

#' @name pep_mac239
#' @title Reference peptide dataset for SIV enveloppe
#'
#' @description A peptide dataset created from the multiple alignment of the 
#' reference sequence mac239 from LANL database, our mac239 and E660 subtypes.
#' 
#' @usage data(pep_mac239)
#' @format A \code{GRanges} object. One row per peptide.
#'
#' @details
#' \itemize{
#'   \item{peptides:}{ The rownames are 15mers from the envelope of mac239 and E660.}
#'   \item{seqnames:}{ An information regarding the genomic location of the peptide.}
#'   \item{ranges:}{ The start and end coordinate of the peptide based on the coordinates of the reference mac239.}
#'   \item{aligned:}{ The aligned column represents the peptide as it is aligned in the multiple alignment.}
#'   \item{peptideNb:}{Number of the corresponding mac239 peptide.}
#'   \item{zSums:}{ The sum of the zScale of each amino-acid of the 15mer.}
#'   \item{clade:}{ A list of the subtypes which contain this peptide.}
#' }
#'
#' @docType data
#' @note The alignments used to create this collections are available in the extdata/alignments folder, as well as the raw sequences on fasta format.
#' @references http://www.hiv.lanl.gov/content/sequence/HIV/REVIEWS/SIV_NUMBERING2001/SivNumbering.html
##
NULL

#' @name pep_hxb2JPT
#' @title Peptide collection for HXB2 and JPT clades
#' 
#' @usage data(pep_hxb2JPT)
#' @format A \code{GRanges} object. One row per peptide.
#'
#' @docType data
##
NULL

#' @name pep_m239smE543
#' @title Peptide collection for SIV enveloppe
#' 
#' @usage data(pep_m239smE543)
#' @format A \code{GRanges} object. One row per peptide.
#'
#' @docType data
##
NULL
