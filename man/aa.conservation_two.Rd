\name{aa.conservation_two}
\alias{aa.conservation_two}
\title{
First and Second Modal Amino Acid Conservation}
\description{
Calculates the conservation of the amino acids that occur first and second most often in a matrix of sequences}
\usage{
aa.conservation_two(modal1, modal2, p, seqlength)
}
\arguments{
  \item{modal1}{
A vector of the frequencies for the amino acids in the first modal sequence.
}
  \item{modal2}{
A vector of the frequencies for the amino acids in the second modal sequence
}
  \item{p}{
A conservation value for the amino acid frequencies to be compared to.
}
  \item{seqlength}{
The length of the amino acid sequence.
}
}
\details{
The argument modal1 can be calculated using the aa.MODE.freq function, and the argument modal2 can be calculated using he aa.MODE.second.freq function.}
\value{
A vector that contains how many amino acids from the first and second modal sequences are conserved at the specified conservation level.
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
aminoAcid_secondModalFreq <- aa.MODE.second.freq(aminoAcid_frequency.Matrix, 216)
aminoAcid_secondConservation_99.9 <- aa.conservation_two(aminoAcid_firstModalFreq, 
    aminoAcid_secondModalFreq, 0.999, 216)}
}
