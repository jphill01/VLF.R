\name{aa.MODE}
\alias{aa.MODE}
\title{
Amino Acid Modal Sequence
}
\description{
Calculates the amino acid sequence that occurs most often in a matrix of amino acid sequences
}
\usage{
aa.MODE(freq.matrix, seqlength)
}
\arguments{
  \item{freq.matrix}{
Freuqncy matrix for amino acids.
}
  \item{seqlength}{
Length of amino acid sequences.
}
}
\details{
The argument freq.matrix can be calculated using the function aa.frequency.matrix.function}
\value{
A vector containing the amino acid sequence that occurs most often.
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{data(birds_aminoAcids)
birds_aminoAcid_count <- aa.count.function(birds_aminoAcids, 216)
aminoAcid_frequency.Matrix <- aa.frequency.matrix.function(birds_aminoAcid_count, 216)
aminoAcid_Modal <- aa.MODE(aminoAcid_frequency.Matrix, 216)}
}