#' Summary tables from pepStat
#' 
#' This is the result tables of a peptide microarray analysis using 
#' \code{pepStat}. It summarizes the antibody binding prediction for each 
#' peptide, depending on the group. \code{restab_aggregate} has one row per 
#' peptide. Peptides that belong to more than one clade have a single entry. 
#' \code{restab} has one row per peptide per clade. Each clade has
#' been normalized separately.
#' 
#' @format A \code{data.frame} containing 1964 rows and 9 variables for 
#' \code{restab}. 1423 rows and 9 variables for \code{restab_aggregate}.
#'  
#' \itemize{
#'   \item \code{peptide}: Peptide sequences
#'   \item \code{position}: The position of peptides on the reference sequence HXB2.
#'   \item \code{space}: The location of the peptide. Here, gp160, the evelope of HIV.
#'   \item \code{start}: The start coordinate of the peptide on the reference sequence.
#'   \item \code{end}: The end coordinate of the peptide on the reference sequence.
#'   \item \code{widt}: The length of the peptides.
#'   \item \code{clade}: The virus subtypes that the peptide belongs to.
#'   \item \code{group1}: Frequency of antibody binding events in the subjects 
#'     of group1 for that peptide.
#'   \item \code{group2}: Frequency of antibody binding events in the subjects
#'     of group2 for that peptide.
#' }
#' 
#' @note For more information, see ?\code{pepStat::restab}.
#' 
#' @name restab
#' @aliases
#' restab
#' restab_aggregate
NULL