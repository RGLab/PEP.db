\docType{data}
\name{pep_hxb2}
\alias{pep_hxb2}
\title{Reference peptide dataset for HIV enveloppe}
\format{A \code{RangedData} object. One row per peptide.}
\usage{
data(pep_hxb2)
}
\description{
A peptide dataset created from the multiple alignment of
the reference sequence hxb2 and the seven subtypes A, B, C,
D, M, CRF01, CRF02
}
\details{
\itemize{ \item{peptides:}{ The rownames are 15mers from
the envelope of 7 subtypes of HIV-1.} \item{space:}{ An
information regarding the genomic location of the peptide.}
\item{ranges:}{ The start and end coordinate of the peptide
based on the coordinates of the reference hxb2.}
\item{aligned:}{ The aligned column represents the peptide
as it is aligned in the multiple alignment.}
\item{trimmed:}{ The trimmed column represents the part of
hxb2's sequence aligned with the peptide.} \item{zSums:}{
The sum of the zScale of each amino-acid of the 15mer.}
\item{clade:}{ A list of the subtypes which contain this
peptide.} }
}
\note{
The alignments used to create this collections are
available in the extdata/alignments folder, as well as the
raw sequences on fasta format.
}
\references{
http://www.hiv.lanl.gov/content/sequence/HIV/REVIEWS/HXB2.html
}

