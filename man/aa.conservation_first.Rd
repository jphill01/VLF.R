\name{aa.conservation_first}
\alias{aa.conservation_first}
\title{
First Modal Amino Acid Conservation
}
\description{
Calculates the conservation of the first amino acids that occur most often in a matrix of amino acid sequences.
}
\usage{
aa.conservation_first(modal, p, seqlength)
}
\arguments{
  \item{modal}{
A vector of the frequencies for the amino acids in the first modal sequence.
}
  \item{p}{
A conservation value for the amino acid frequencies to be compared to
}
  \item{seqlength}{
Length of the amino acid sequence
}
}
\details{
The item modal can be calculated using the aa.MODE.freq function.
}
\value{
A vector that contains how many amino acids from the first modal sequence are conserved at the specified conservation level.
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{data(birds_aminoAcids)
birds_aminoAcid_count <- aa.count.function(birds_aminoAcids, 216)
aminoAcid_frequency.Matrix <- aa.frequency.matrix.function(birds_aminoAcid_count, 216)
aminoAcid_Modal <- aa.MODE(aminoAcid_frequency.Matrix, 216)
aminoAcid_firstModalFreq <- aa.MODE.freq(aminoAcid_frequency.Matrix, 216)
aminoAcid_firstConservation_100 <- aa.conservation_first(aminoAcid_firstModalFreq, 1, 216)}
}