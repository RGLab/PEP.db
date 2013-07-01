\name{pep_m239smE543}
\alias{pep_m239smE543}

\title{
	Peptide dataset for SIV envelope.
}

\description{
A peptide dataset created from the multiple alignment of the reference sequence mac239 from LANL database, our mac239 and E43 subtypes.

\describe{
\item{peptides:}{The rownames are 15mers from the envelope of mac239 and E543.}
\item{space:}{An information regarding the genomic location of the peptide.}
\item{ranges:}{The start and end coordinate of the peptide based on the coordinates of the reference hxb2.}
\item{aligned:}{The aligned column represents the peptide as it is aligned in the multiple alignment.}
\item{zSums:}{The sum of the zScale of each amino-acid of the 15mer.}
\item{clade:}{A list of the subtypes which contain this peptide.}
}
}
\note{
	The alignments used to create this collections are available in the extdata/alignments folder, as well as the raw sequences on fasta format.
}

\keyword{datasets}
